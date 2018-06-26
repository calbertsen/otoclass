##' @importFrom Matrix sparseMatrix diag
contourQ <- function(N, delta = 0){
    pnts <- c(1,1,2)
    t <- seq(0,1,len = N)
    for(i in 2:length(t)){
        pnts <- rbind(pnts,
                      c(i,i,2),
                      c(i,i-1,-1),
                      c(i-1,i,-1)
                      )
    }
    pnts <- rbind(pnts,
                  c(length(t),1,-1),
                  c(1,length(t),-1)
                  )
    Q <- Matrix::sparseMatrix(i=pnts[,1],j=pnts[,2],x=pnts[,3])
    Matrix::diag(Q) <- Matrix::diag(Q) + delta
    return(Q)
}


##' Maximum Likelihood Contour Discrimination
##'
##' @param train 3D Array of training data (2 x number of contour points x number of observations) 
##' @param group Factor of training groups
##' @param test 3D array of test data
##' @param prior Prior probability of groups
##' @param penalty p to use for Lp penalty. Zero is no penalty.
##' @param lambda Positive scalar factor for Lp penalty. Zero is no penalty.
##' @param Nefd Number of EFDs to use for mean otolith
##' @param correlatedCoordinates Should coordinates be correlated?
##' @param equalVariance Should the variance be equal in the two coordinates?
##' @param silent Should the TMB object be silent?
##' @param control control parameters passes to nlminb
##' @param ... Other parameters
##' @return a list of the result
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats runif nlminb
##' @importFrom TMB MakeADFun
##' @export
mlcd <- function(train, group, test,
                 prior = as.vector(table(group)) / length(group),
                 penalty = 0,
                 lambda = 0.4,
                 Nefd = 60,
                 correlatedCoordinates = FALSE,
                 equalVariance = TRUE,
                 silent = FALSE,                 
                 control = list(iter.max = 100000, eval.max = 100000),
                 ...){
##### Checks #####
    ## 1) test/train should be an array with three dimensions
    ## 2) first dimension of test/train should be 2
    ## 3) Length of group should be same as number of test samples (dim 3)
    ## 4) Group should be a factor
    if(!is.factor(group))
        group <- factor(group)
    ## 5) prior should be same length as number of groups

##### Prepare Q #####
    N <- dim(train)[2]                 
    Q <- contourQ(N,0)


##### Data for TMB #####
    dat <- list(model = 2L,
                X = train,
                G = group,
                Q = Q,
                penalty = penalty,
                prior = prior,
                X_pred = test,
                logLambda = log(lambda)
                )
##### Parameters for TMB #####
    ## Efds
    n <- nrow(dat$X)
    efds <- array(0,dim = c(4,Nefd,nlevels(dat$G)))
    efds[] <- stats::runif(prod(dim(efds)),-0.01,0.01)
    efds[1,1,] <- 1
    efds[2:3,1,] <- 0
    efds[4,1,] <- 1
    ## par list
    par <- list(mu = matrix(0,0,0),
                efd = efds,
                logSigma = matrix(0,n,nlevels(dat$G)),
                corpar = matrix(0,(n*n-n)/2,nlevels(dat$G)),
                logDelta = numeric(nlevels(dat$G))
                )
##### Map for TMB #####
    ## normalized efd: in first harmonic 3 variables are degenerate
    efdMap <- array(1:prod(dim(efds)),dim = dim(efds))
    efdMap[1:3,1,] <- NA
    ## Corpar map
    if(correlatedCoordinates){
        corparMap <- factor(1:length(par$corpar))
    }else{
        corparMap <- factor(rep(NA,length(par$corpar)))
    }
    ## variance parameter map
    if(equalVariance){
        sigmaMap <- factor(row(par$logSigma))
    }else{
        sigmaMap <- factor(1:length(par$logSigma))
    }
    map <- list(logSigma = sigmaMap,
                corpar = corparMap,
                efd = factor(efdMap)
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
