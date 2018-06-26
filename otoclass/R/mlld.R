##' Maximum Likelihood Linear Discrimination
##'
##' Equivalent to LDA
##' @param train Matrix of training data (number of observations x number of features)
##' @param group Factor of training groyps
##' @param test Matrix of test data
##' @param prior Prior probability of groups
##' @param penalty p to use for Lp penalty. Zero is no penalty
##' @param lambda Positive scalar factor for Lp penalty. Zero is no penalty.
##' @param independent Should features be treated as independent?
##' @param silent Should the TMB object be silent?
##' @param control control parameters passes to nlminb
##' @param ... Other parameters
##' @return a list of the result
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats nlminb cor
##' @importFrom methods as
##' @importFrom TMB MakeADFun
##' @export
mlld <- function(train, group, test,
                 prior = as.vector(table(group)) / length(group),
                 penalty = 0,
                 lambda = 0.4,
                 independent = FALSE,
                 silent = FALSE,                 
                 control = list(iter.max = 100000, eval.max = 100000),
                 ...){
##### Checks #####
    ## 1) test/train should be a matrix
    ## 2) first dimension of test/train should be same as length of group)
    ## 3) Group should be a factor
    if(!is.factor(group))
        group <- factor(group)
    ## 4) prior should be same length as number of groups

##### Prepare Q #####
    N <- dim(train)[2]                 
    Q <- methods::as(matrix(0.0,0,0),"sparseMatrix")

##### Data for TMB #####
    dat <- list(model = ifelse(independent,0L,1L),
                X = t(train),
                G = group,
                Q = Q,
                penalty = penalty,
                prior = prior,
                X_pred = t(test),
                logLambda = log(lambda)
                )
##### Parameters for TMB #####
    n <- nrow(dat$X)
    ## par list
    mn <- sapply(levels(dat$G),function(i)apply(dat$X[,dat$G==i,drop=FALSE],1,mean))
    if(independent){
        corcalc <- 0
    }else{
        corcalc <- t(chol(stats::cor(t(dat$X))))[lower.tri(t(chol(stats::cor(t(dat$X)))),diag = FALSE)]
    }
    par <- list(mu = mn, ##matrix(0.0,nrow(dat$X),nlevels(dat$G)),
                efd = matrix(0.0,0,0),
                logSigma = matrix(log(apply(dat$X,1,sd))+2,n,nlevels(dat$G)),
                corpar = matrix(corcalc,(n*n-n)/2,nlevels(dat$G)),
                logDelta = 0
                )
##### Map for TMB #####
    ## normalized efd: in first harmonic 3 variables are degenerate
    ## variance parameter map
    if(independent){
        corparMap <- factor(rep(NA,length(par$corpar)))
    }else{
        corparMap <- factor(row(par$corpar))
    }
    sigmaMap <- factor(row(par$logSigma))
    map <- list(logSigma = sigmaMap,
                corpar = corparMap,
                logDelta = factor(NA)
                )

##### Make TMB object #####
    obj <- TMB::MakeADFun(dat,par,map,
                          silent = silent,
                          DLL = "otoclass")
    opt <- stats::nlminb(obj$par,
                         obj$fn,
                         obj$gr,
                         control = control)
    rp <- obj$report(opt$par)
    

    return(list(obj = obj, opt = opt, rp = rp))
    
}
