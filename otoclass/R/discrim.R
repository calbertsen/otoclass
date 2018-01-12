
mvnorm <- function(x,mn,sig){
    -0.5*log(det(sig)) - 0.5*(x-mn)%*%solve(sig)%*%(x-mn) -0.5*dim(sig)[1]*log(2*pi)
}
which.max.safe <- function(x,grplevel){
    if(all(is.na(x)))
        return(NA)
    grplevel[which.max(x)]
    
}

#' Discriminant analysis
#'
#' \code{discrim} performs linear and quadratic discriminant analysis.
#'
#' @param train A matrix of (continuous) features for the training sample
#' @param group A vector of group labels for the training sample
#' @param test A matrix of (continuous) features for the test sample
#' @param type Type of discriminant analysis to use. lda for linear and qda for quadratic.
#' @param verbose Write maximum gradient components to the terminal?
#' @param dist Distribution to use (currently ignored)
#'
#' @return stuff
#' @importFrom stats cov
#' @export
discrim <- function(train, group, test,type = "lda", verbose=TRUE, dist = "normal", prior = table(group)/length(group), nrr = rank){

    
    if(type == "rrlda"){
        mn <- sapply(1:nlevels(group),
                     function(i)apply(train[as.numeric(group)==i,],2,mean))
        W <- stats::cov(train-t(mn[,group]))
        B <- stats::cov(apply(t(mn)-apply(train,2,mean),2,function(x)x*prior))
        ew <- svd(W)
        W_sqrt <- diag(sqrt(ew$d))%*%t(ew$v)
        W_msqrt <- solve(W_sqrt)
        Bstar <- t(W_msqrt)%*%B%*%W_msqrt
        ebs <- eigen(Bstar)
        rank <- sum(ebs$values>1e-04)
        transfmat <- W_msqrt%*%ebs$vectors[,1:nrr]
        testUse <- test%*%transfmat
        trainUse <- train%*%transfmat

    }else{
        transfmat <- diag(1,dim(train)[2])
        testUse <- test
        trainUse <- train
    }

    lpi <- log(prior)

 
    mn <- sapply(1:nlevels(group),
                 function(i)apply(matrix(trainUse[as.numeric(group)==i,],ncol=dim(trainUse)[2]),2,mean))

    mn <- matrix(mn,ncol=nlevels(group))
    
    if(type == "qda"){
        sig <- sapply(1:nlevels(group),
                      function(i)list(cov(matrix(trainUse[as.numeric(group)==i,]-mn[,i],ncol=dim(trainUse)[2]))))

        lnorm <- function(i){
            apply(testUse,1,function(x)mvnorm(as.vector(x),mn=as.vector(mn[,i]),sig=sig[[i]]))+lpi[i]
        }

    }else{
        sig <- cov(trainUse-t(matrix(mn[,group],nrow=dim(mn)[1])))
        lnorm <- function(i){
            apply(testUse,1,function(x)mvnorm(as.vector(x),mn=as.vector(mn[,i]),sig=sig))+lpi[i]
        }
    }

    lps <- sapply(1:nlevels(group),lnorm)
    if(nrow(testUse) == 1)
        lps <- matrix(lps,1)
    prop <- t(apply(lps,1,function(x)exp(x)/sum(exp(x))))
    prop[is.nan(prop)] <- NA
    res <- list()
    class(res) <- "oto_discrim"

    res$type <- type
    res$probabilities <- prop
    colnames(res$probabilities) <- levels(group)
    #rownames(res$probabilities) <- rownames(test)
    #res$predicted <- factor(levels(group)[apply(prop,1,which.max)], levels = levels(group))
    res$predicted <- factor(apply(prop,1,which.max.safe,grplevel=levels(group)),levels = levels(group))
 
    res$tranfmat <- transfmat
    res$sig <- sig
    if(type=="qda"){
        res$sig <- lapply(res$sig,function(x){y <- x;colnames(y) <- rownames(y) <- colnames(train);y})
        names(res$sig) <- levels(group)
    }
    colnames(mn) <- levels(group)
    rownames(mn) <- colnames(train)
    res$mn <- mn
    return(res)
}











##' discrim with TMB
##'
##' @param train 
##' @param group 
##' @param test 
##' @param type 
##' @param verbose 
##' @param dist 
##' @param prior 
##' @return ...
##' @author Christoffer Moesgaard Albertsen
##' @importFrom TMB MakeADFun
##' @importFrom stats nlminb
discrimTMB <- function(train, group, test,type = "lda", verbose=TRUE, dist = "normal", prior = table(group)/length(group)){

    if(!is.factor(group)){
        group <- factor(group)
    }
    
    dat <- list(group = group,
                covariates = train,
                pi = prior,
                distribution = 1,
                predCovariates = test)
    param <- list(means = matrix(0,
                      nrow = dim(dat$covariates)[2],
                      ncol = nlevels(dat$group)),
                  logSd = matrix(0,
                      nrow = dim(dat$covariates)[2],
                      ncol = nlevels(dat$group)),
                  trho  = matrix(0,
                      nrow = 0.5*dim(dat$covariates)[2]*(dim(dat$covariates)[2]-1),
                      ncol = nlevels(dat$group))
                  )

    if(type=="lda"){
        map <- lda_map(param)
    }else if(type=="qda"){
        map <- list()
    }else{
        stop("Invalid type. Must be lda or qda.")
    }


    
    obj <- TMB::MakeADFun(dat,param,map=map,DLL="discrim",silent=!verbose)
    opt <- stats::nlminb(obj$par,obj$fn,obj$gr,obj$he)
    
    res <- list()
    class(res) <- "oto_discrim"

    res$type <- type
    res$distribution <- "normal"
    res$probabilities <- obj$report()$probPred
    colnames(res$probabilities) <- levels(group)
    rownames(res$probabilities) <- rownames(test)
    res$predicted <- factor(levels(group)[apply(obj$report()$probPred,1,which.max)], levels = levels(group))
    res$logLik <- opt$value
    res$convergence <- opt$convergence
    res$means <- obj$env$parList(opt$par)$means
    colnames(res$means) <- levels(group)
    rownames(res$means) <- colnames(test)
    res$covariance <- obj$report(opt$par)$covs
    dimnames(res$covariance) <- list(colnames(test),colnames(test),levels(group))

    return(res)
    
}
