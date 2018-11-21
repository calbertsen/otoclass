##' Maximum Likelihood Linear Discrimination
##'
##' Equivalent to LDA
##' @param train Matrix of training data (number of observations x number of features)
##' @param group Factor of training groyps
##' @param test Matrix of test data
##' @param prior Prior probability of groups
##' @param penalty p to use for Lp penalty. Zero is no penalty
##' @param lambda Positive numeric vector factor for Lp penalty. Zero is no penalty. The first element determines the penalty on group differences, the second element determines the penalty on absolute parameter values, and the third element determines the penalty on correlation parameters.
##' @param independent Should features be treated as independent?
##' @param silent Should the TMB object be silent?
##' @param control control parameters passes to nlminb
##' @param formula Formula to account for in the linear discrimination
##' @param data Covariates for training data
##' @param dataTest Covariates for testing data
##' @param ... Other parameters
##' @return a list of the result
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats nlminb cor
##' @importFrom methods as
##' @importFrom TMB MakeADFun sdreport
##' @export
mlld <- function(train, group, test,
                 prior = as.vector(table(group)) / length(group),
                 featureSelection = c("none","fcbf","Lp"),
                 penalty = 0,
                 lambda = 0.4,
                 estimateLambda = TRUE,
                 independent = FALSE,
                 silent = FALSE,                 
                 control = list(iter.max = 100000, eval.max = 100000),
                 formula = ~1,
                 formulaCommon = ~1,
                 data = NULL,
                 dataTest = NULL,
                 estimateUnbiasedTestProportions = FALSE,
                 trainContributionToTestProportions = FALSE,
                 biasCorrectionGroup = NULL,
                 trainBiasCorrectionGroup = NULL,
                 addMisfitCategory = FALSE,
                 fcbfDelta = 0.01,
                 fcbfLambda = 1,
                 onlyObj = FALSE,
                 doSdreport = FALSE,
                 ...){
##### Checks #####
    featureSelection <- match.arg(featureSelection)
    call <- match.call()
    ##args <- c(as.list(environment()), list(...))
    if(missing(test))
        test <- train
    ## 1) test/train should be a matrix
    ## 2) second dimension of test/train should be same as length of group)
    ## 3) Group should be a factor
    if(!is.factor(group))
        group <- factor(group)
    ## 4) prior should be same length as number of groups
    ## if(length(lambda) != 3){
    ##     warning(sprintf("lambda should be a vector of length 3. %s",
    ##                     ifelse(length(lambda) > 3,
    ##                            "The values have been repeated.",
    ##                            "The first three elements are used.")
    ##                     ))
    ##     lambda <- rep(lambda,length = 3)
    ## }        
    if(length(lambda) > 1 | lambda <= 0)
        stop("Lambda must be a positive scalar.")
    if(!identical(formula, ~1) & (is.null(data) | is.null(dataTest))){
        stop("When a formula is specified, data and dataTest must be given.")
    }
    
##### Prepare Q #####
    N <- dim(train)[2]                 
    Q <- methods::as(matrix(0.0,0,0),"sparseMatrix")

##### Prepare covar #####
    if(is.null(data))
        data <- data.frame(ID = 1:nrow(train))
    mf <- model.frame(formula, data)
    
    if(is.null(dataTest))
            dataTest <- data.frame(ID = 1:nrow(test))
    if(!identical(nrow(train), nrow(data)))
        stop("data and train must have the same number of rows.")
    if(!identical(nrow(test), nrow(dataTest)))
        stop("dataTest and test must have the same number of rows.")
    covar <- Matrix::sparse.model.matrix(terms(mf),
                                         data = data,
                                         transpose = TRUE,
                                         row.names = FALSE)
    if(inherits(covar,"dgCMatrix"))
        covar <- as(covar,"dgTMatrix")
    covarTest <- Matrix::sparse.model.matrix(terms(mf),
                                             data = dataTest,
                                             transpose = TRUE,
                                             row.names = FALSE)
    if(inherits(covarTest,"dgCMatrix"))
        covarTest <- as(covarTest,"dgTMatrix")
    indxUse <- switch(featureSelection,
                      none = 1:ncol(train),
                      fcbf = fcbf(train, group, delta = fcbfDelta, lambda = fcbfLambda),
                      Lp = 1:ncol(train),
                      default = 1:ncol(train)
                      )

    if(is.null(biasCorrectionGroup)){
        biasCorrectionGroup <- factor(rep(1,nrow(test)))
    }else{
        if(!(is.factor(biasCorrectionGroup) & length(biasCorrectionGroup) == nrow(test)))
            stop("biasCorrectionGroup must be a factor of the same length as test")
    }

    ## Common effects
    formulaCommon <- update.formula(formulaCommon, ~ . +1)
    mfCommon <- model.frame(formulaCommon, data)
    commonCovar <- Matrix::sparse.model.matrix(terms(mfCommon),
                                         data = data,
                                         transpose = TRUE,
                                         row.names = FALSE)
    if(inherits(commonCovar,"dgCMatrix"))
        commonCovar <- as(commonCovar,"dgTMatrix")
   commonCovarTest <- Matrix::sparse.model.matrix(terms(mfCommon),
                                             data = dataTest,
                                             transpose = TRUE,
                                             row.names = FALSE)
    if(inherits(commonCovarTest,"dgCMatrix"))
        commonCovarTest <- as(commonCovarTest,"dgTMatrix")

    if(!is.null(trainBiasCorrectionGroup)){
        if(trainContributionToTestProportions){
            if(!all(levels(biasCorrectionGroup) == levels(trainBiasCorrectionGroup)))
                stop("biasCorrectionGroup and trainBiasCorrectionGroup must have the same levels")
            if(length(trainBiasCorrectionGroup) != nrow(train))
                stop("trainBiasCorrectionGroup must have the same length as the number of rows in train")
        }
    }else{      
        trainBiasCorrectionGroup <- factor(rep(1,nrow(train)))
    }
    
##### Data for TMB #####
    dat <- list(model = ifelse(independent,0L,1L),
                X = t(train[,indxUse, drop = FALSE]),
                G = group,
                Q = Q,
                penalty = penalty * (featureSelection == "Lp"),
                prior = prior,
                X_pred = t(test[,indxUse, drop = FALSE]),
                covar = covar,
                covar_pred = covarTest,
                commonCovar = commonCovar,
                commonCovar_pred = commonCovarTest,
                estimateUnbiasedTestProportions = as.integer(estimateUnbiasedTestProportions),
                trainContributionToTestProportions = as.integer(trainContributionToTestProportions),               testBiasCorrectionGroup = biasCorrectionGroup,
                trainBiasCorrectionGroup = trainBiasCorrectionGroup,
                addMisfitCategory = as.integer(addMisfitCategory)
                )
##### Parameters for TMB #####
    n <- nrow(dat$X)
    ## par list
    mndim <- c(nrow(dat$covar), nrow(dat$X), nlevels(dat$G))
    mn <- array(rnorm(prod(mndim),0,0.001), dim = mndim)
    ##mn <- sapply(levels(dat$G),function(i)apply(dat$X[,dat$G==i,drop=FALSE],1,mean))
    if(independent){
        corcalc <- 0
    }else{
        corcalc <- t(chol(stats::cor(t(dat$X))))[lower.tri(t(chol(stats::cor(t(dat$X)))),diag = FALSE)]
    }
    if(!is.null(call$parlist)){
        par <- eval(call$parlist)
        if(!is.null(call$lambda))
            par$logLambda <- log(lambda)
    }else{
        par <- list(mu = mn, ##matrix(0.0,nrow(dat$X),nlevels(dat$G)),
                    commonMu = matrix(0,nrow(dat$commonCovar),nrow(dat$X)),
                    efd = matrix(0.0,0,0),
                    logSigma = matrix(log(apply(dat$X,1,sd))+2,n,nlevels(dat$G)),
                    corpar = matrix(corcalc,(n*n-n)/2,nlevels(dat$G)),
                    logDelta = 0,
                    logLambda = log(lambda),
                    thetaIn = matrix(0.0,(nlevels(dat$G)-1 + as.numeric(addMisfitCategory))*estimateUnbiasedTestProportions,
                                     nlevels(biasCorrectionGroup)*estimateUnbiasedTestProportions)
                    )
    }
    np <- Reduce("+",lapply(par,length))
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

    if(nrow(par$commonMu) > 0){
        cMuMap <- matrix(1:length(par$commonMu),nrow(par$commonMu),ncol(par$commonMu))
        cMuMap[1,] <- NA    
        map$commonMu <- factor(cMuMap)
    }
    if(!estimateLambda || featureSelection != "Lp"){
        map$logLambda <- factor(par$logLambda*NA)
        rnd <- c()
    }else{
        rnd <- c("mu","commonMu")
    }
    npfix <- Reduce("+",lapply(map,function(x)sum(is.na(x))))
    if(np-npfix > nrow(train) / 10)
        warning(sprintf("The model has %s parameters with only %s observations - %.2f observations per parameter. There is a high risk of overfitting.",np-npfix,nrow(train),nrow(train)/(np-npfix))) 
    
##### Make TMB object #####
    obj <- TMB::MakeADFun(dat,par,map,
                          silent = silent, profile = rnd,
                          DLL = "otoclass")
    if(onlyObj)
        return(obj)
    opt <- stats::nlminb(obj$par,
                         obj$fn,
                         obj$gr,
                         control = control)
    rp <- obj$report(obj$env$last.par.best)

    muNames <- list(rownames(dat$covar),
                    colnames(train),
                    levels(group))
    commonMuNames <- list(rownames(dat$commonCovar),
                          colnames(train))

    posterior <- exp(t(rp$logpred))
    colnames(posterior) <- levels(group)
    rownames(posterior) <- rownames(test)
    clas <- factor(levels(group)[apply(posterior,1,which.max)], levels = levels(group))

    if(doSdreport){
        sdr <- TMB::sdreport(obj,opt$par)
    }else{
        sdr <- NULL
    }
    
    res <- list(tmb_data = obj$env$data,
                tmb_map = map,
                tmb_random = rnd,
                X_cols = ncol(train),
                all_vars = colnames(get_all_vars(terms(mf), data)),
                silent = silent,
                opt = opt,
                rp = rp,
                sdr = sdr,
                featureSelection = featureSelection,
                features = indxUse,
                overfit = np-npfix > nrow(train) / 10,
                muNames = muNames,
                commonMuNames = commonMuNames,
                pl = obj$env$parList(opt$par),
                call = call,
                prior = prior,
                terms = terms(mf),
                termsCommon = terms(mfCommon),
                train = train,
                test = test,
                group = group,
                data = data,
                test_class = clas,
                test_posterior = posterior,
                varBetweenGroups = NA,
                varWithinGroups = NA
                )
    class(res) <- "mlld"
    res$varBetweenGroups <- calculateVarBetweenGroups(res)
    res$varWithinGroups <- calculateVarWithinGroups(res)
    
    return(res)    
}

##' @export
print.mlld <- function(x, ...){
    print(x$opt)
}
    
