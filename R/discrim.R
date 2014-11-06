
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


#' @export
discrim <- function(train, group, test,type = "lda", verbose=TRUE, dist = "normal", prior = table(group)/length(group), nrr = rank){

    lpi <- log(prior)

    if(type == "rrlda"){
        G <-  model.matrix(~group-1)
        M <- t(sapply(1:nlevels(group),
                      function(i)apply(train[as.numeric(group)==i,],2,mean)))
        W <- t(train-G%*%M)%*%(train-G%*%M)/(-diff(dim(G)))
        ew <- eigen(W)
        Dsqinv <- diag(sqrt(1/ew$values))
        trainS <- train%*%ew$vectors%*%Dsqinv
        MS <- t(sapply(1:nlevels(group),
                       function(i)apply(trainS[as.numeric(group)==i,],2,mean)))
        WS <- t(trainS-G%*%MS)%*%(trainS-G%*%MS)/(-diff(dim(G)))

        xbar <- prior%*%t(MS)
        One <- rep(1,dim(G)[1])
        Ct <- t(trainS-One%*%xbar)%*%(trainS-One%*%xbar)/(dim(G)[1]-dim(xbar)[2])
        B <- Ct-WS
        eb <- eigen(B)
        rank <- sum(eb$values/sum(eb$values)>1e-03)
        transfmat <- ew$vectors%*%Dsqinv%*%eb$vectors[,1:min(nrr,
                                                             rank,
                                                             nlevels(group)-1)]
        propExpl <- (eb$values/sum(eb$values))[1:min(nrr,
                                                     rank,
                                                     nlevels(group)-1)]
        testUse <- test%*%transfmat
        trainUse <- train%*%transfmat

        mn <- t(M%*%transfmat)
        sig <- diag(rep(1,dim(mn)[1]))
        lnorm <- function(i){
            apply(testUse,1,function(x)mvnorm(as.vector(x),
                                              mn=as.vector(mn[,i]),
                                              sig=sig))+lpi[i]
        }
        lps <- sapply(1:nlevels(group),lnorm)
    }else{
        transfmat <- diag(1,dim(train)[2])
        testUse <- test
        trainUse <- train
        propExpl <- NULL

       
 
        mn <- sapply(1:nlevels(group),
                     function(i)apply(matrix(trainUse[as.numeric(group)==i,],
                                             ncol=dim(trainUse)[2]),
                                      2,mean))

        mn <- matrix(mn,ncol=nlevels(group))
    
        if(type == "qda"){
            sig <- sapply(1:nlevels(group),
                          function(i)list(
                              cov(matrix(trainUse[as.numeric(group)==i,]-mn[,i],
                                         ncol=dim(trainUse)[2]))))

            lnorm <- function(i){
                apply(testUse,1,function(x)mvnorm(as.vector(x),
                                                  mn=as.vector(mn[,i]),
                                                  sig=sig[[i]]))+lpi[i]
            }

        }else{
            sig <- cov(trainUse-t(matrix(mn[,group],nrow=dim(mn)[1])))
            lnorm <- function(i){
                apply(testUse,1,function(x)mvnorm(as.vector(x),
                                                  mn=as.vector(mn[,i]),
                                              sig=sig))+lpi[i]
            }
        }

        lps <- sapply(1:nlevels(group),lnorm)
    }
    prop <- t(apply(lps,1,function(x)exp(x)/sum(exp(x))))
    prop[is.nan(prop)] <- NA
    res <- list()
    class(res) <- "oto_discrim"

    res$type <- type
    res$logProbabilities <- lps
    res$probabilities <- prop
    colnames(res$probabilities) <- levels(group)
    #rownames(res$probabilities) <- rownames(test)
    #res$predicted <- factor(levels(group)[apply(prop,1,which.max)], levels = levels(group))
    res$predicted <- factor(apply(lps,1,which.max.safe,grplevel=levels(group)),levels = levels(group))
 
    res$tranfmat <- transfmat
    res$sig <- sig
    res$mn <- mn
    res$proportionExplained <- propExpl
    return(res)
}












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


    
    obj <- MakeADFun(dat,param,map=map,DLL="discrim")
    obj$env$inner.control$trace <- verbose
    obj$env$tracemgc <- verbose

    opt <- nlminb(obj$par,obj$fn,obj$gr,obj$he)
    
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
