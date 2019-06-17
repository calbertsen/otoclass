##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param proportionGroup 
##' @param confusionGroup 
##' @param ... 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats nlminb cor model.frame terms update.formula qlogis .getXlevels get_all_vars
##' @importFrom methods as
##' @importFrom TMB MakeADFun sdreport
##' @importFrom expm logm
##' @export
mlld <- function(## Data related
                 y,
                 group,
                 data = NULL,
                 formula = ~1,
                 formulaCommon = ~1,
                 proportionGroup = factor(rep("Baseline",nrow(y))),
                 confusionGroup = factor(ifelse(is.na(group),"Unknown","Known")), ##factor(rep("Known",nrow(y))),
                 ## Penalty related
                 lp_penalty = 0,
                 lambda = 0.4,
                 estimateLambda = TRUE,
                 ## TMixture related
                 tMixture = 0,
                 tDf = 5,
                 estimateTMix = FALSE,
                 sameTMix = FALSE,
                 estimateTDf = FALSE,
                 sameTDf = FALSE,
                 ## KW distribution related
                 aParam = 1,
                 bParam = 1,
                 estimateKWAParam = FALSE,
                 estimateKWBParam = FALSE,
                 useSecondKWParameterization = FALSE, ##estimateKWAParam && estimateKWBParam,
                 ## Settings
                 independent = FALSE,
                 silent = FALSE,                 
                 control = list(iter.max = 100000, eval.max = 100000),
                 drop.unused.levels = TRUE,
                 onlyObj = FALSE,
                 doSdreport = FALSE,
                 getReportCovariance = FALSE,
                 equalVariance = TRUE,
                 confusionMatrixArray = NULL,
                 confusionLevelTypes = levels(confusionGroup),  ## Possible values: c("Known","Fixed","Estimate","Unknown")               
                 ...){
    
    cl <- match.call()
##### Checks #####
    if(!is.matrix(y))
        if(is.vector(y)){
            y <- matrix(y,ncol = 1)
        }else{
            stop("y must be a matrix of observations")
        }
    if(!is.factor(group))
        group <- factor(group)
    if(!(length(lambda) == 2 || length(lambda)==1) | any(lambda <= 0))
        stop("Lambda must be a positive scalar.")
    ## if((!identical(formula, ~1) | !identical(formulaCommon, ~1)) & is.null(data))
    ##     stop("When a formula is specified, data must be given.")
    if(!is.null(data) & !identical(nrow(y), nrow(data)))
        stop("data and y must have the same number of rows.")

    confusionLevelTypes <- c("Known","Fixed","Estimate","Unknown")[pmatch(confusionLevelTypes,c("Known","Fixed","Estimate","Unknown"), duplicates.ok = TRUE)]
    if(any(is.na(confusionLevelTypes)))
        stop('confusionLevelTypes must be a vector of "Known","Fixed","Estimate","Unknown". Partial matches are allowed.')
    if((length(confusionLevelTypes) != nlevels(confusionGroup)))
        stop("The length of confusionLevelTypes must match the number of levels in confusionGroup or the number of fixed levels.")
    if(!is.null(confusionMatrixArray) && !((dim(confusionMatrixArray)[3] == sum(confusionLevelTypes=="Fixed")) || (dim(confusionMatrixArray)[3] == length(confusionLevelTypes))))
        stop("...")

##### Prepare confusion matrix array #####
    hasAll <- !is.null(confusionMatrixArray) && (dim(confusionMatrixArray)[3] == length(confusionLevelTypes))
    CMA <- array(NA,
                 dim = c(nlevels(group), nlevels(group), nlevels(confusionGroup)))
    if(hasAll){
        CMA <- confusionMatrixArray
    }else{
        for(i in 1:nlevels(confusionGroup)){
            if(confusionLevelTypes[i] == "Known"){
                CMA[,,i] <-log(diag(1,nrow = dim(CMA)[1]))
            }else if(confusionLevelTypes[i] == "Fixed"){
                isFixed <- confusionLevelTypes[1:i]=="Fixed"
                CMA[,,i] <- log(confusionMatrixArray[,,sum(isFixed)])
            }else if(confusionLevelTypes[i] == "Estimate"){
                CMA[,,i] <- log(matrix(0.1/(dim(CMA)[1]-1),nrow = dim(CMA)[1], ncol = dim(CMA)[2]))
                diag(CMA[,,i]) <- log(0.9)
            }else if(confusionLevelTypes[i] == "Unknown"){
                CMA[,,i] <- log(matrix(1,nrow = dim(CMA)[1], ncol = dim(CMA)[2]))
            }
        }
    }

    CMA_map <- array(1:prod(dim(CMA)), dim = dim(CMA))
    for(i in 1:dim(CMA_map)[3]){
        if(confusionLevelTypes[i] %in% c("Estimate")){
            diag(CMA_map[,,i]) <- NA
        }else{
            CMA_map[,,i] <- NA
        }
    }
    
  
##### Prepare model matrix #####
    if(is.null(data) & identical(formula,~1) & identical(formulaCommon, ~1))
        data <- data.frame(ID = 1:nrow(y))
    if(!is.null(data) & !is.data.frame(data))
        data <- as.data.frame(data)

    mf <- stats::model.frame(formula, data)

    X <- Matrix::sparse.model.matrix(stats::terms(mf),
                                     data = mf,
                                     transpose = TRUE,
                                     row.names = FALSE,
                                     drop.unused.levels = drop.unused.levels)
    if(inherits(X,"dgCMatrix"))
        X <- as(X,"dgTMatrix")

##### Prepare common model matrix #####
    formulaCommon <- stats::update.formula(formulaCommon, ~ . +1)
    mfCommon <- stats::model.frame(formulaCommon, data)
    XCom <- Matrix::sparse.model.matrix(stats::terms(mfCommon),
                                        data = mfCommon,
                                        transpose = TRUE,
                                        row.names = FALSE,
                                        drop.unused.levels = drop.unused.levels)
    if(inherits(XCom,"dgCMatrix"))
        XCom <- as(XCom,"dgTMatrix")

##### Prepare confusion matrices #####
    MIn <- CMA
    ## for(i in 1:dim(MIn)[3]){
    ##     if(dim(MIn)[1] == 1){
    ##         MIn[,,i] <- 1
    ##     }else{
    ##         MIn[,,i] <- expm::logm(MIn[,,i])
    ##         diag(MIn[,,i]) <- 1
    ##     }
    ##     MIn[,,i] <- log(MIn[,,i])
    ## }

##### Handle NA in group #####
    Guse <- group
    if(any(is.na(Guse))){
        indx <- which(is.na(Guse))
        if(any(confusionLevelTypes[confusionGroup[indx]] != "Unknown"))
            stop("All observations with NA in group must have a confusionGroup with confusionLevelType equal to 'Unknown'")
        Guse[indx] <- levels(group)[1]
    }
    
##### Data for TMB #####
    dat <- list(model = ifelse(independent,0L,1L),
                G = Guse,
                Y = t(y),
                X = X,
                XCom = XCom,
                proportionGroup = proportionGroup,
                confusionGroup = confusionGroup,
                penalty = lp_penalty,
                Y_pred = matrix(0,0,0),
                X_pred = as(matrix(0,0,0),"dgTMatrix"),
                XCom_pred = as(matrix(0,0,0),"dgTMatrix"),
                proportionGroup_pred = factor(levels=levels(Guse)),
                kw_idparam = as.integer(useSecondKWParameterization)
                )
    
##### Parameters for TMB #####
    n <- nrow(dat$Y)
    ## par list
    mndim <- c(nrow(dat$X), nrow(dat$Y), nlevels(dat$G))
    mn <- array(0.1, dim = mndim)
    ##mn <- sapply(levels(dat$G),function(i)apply(dat$X[,dat$G==i,drop=FALSE],1,mean))
    if(independent){
        corcalc <- 0
    }else{
        ##corcalc <- t(chol(stats::cor(t(dat$Y))))[lower.tri(t(chol(stats::cor(t(dat$Y)))),diag = FALSE)]
        corcalc <- 0
    }
    if(!is.null(cl$parlist)){
        par <- eval(cl$parlist)
        if(!is.null(cl$lambda))
            par$logLambda <- log(lambda)
    }else{
        sd2 <- function(x){
            if(length(x) == 1)
                return(1)
            return(sd(x))
        }
        if(estimateTMix && missing(tMixture)){
            tMixture <- 0.05
        }
        if(estimateTMix){
            tMixture <- (1.0 - .Machine$double.eps) * (tMixture - .5) + .5
        }
        par <- list(mu = mn,
                    commonMu = matrix(0,nrow(dat$XCom),nrow(dat$Y)),
                    logSigma = matrix(log(apply(dat$Y,1,sd2))+2,n,nlevels(dat$G)),
                    corpar = matrix(corcalc,(n*n-n)/2,nlevels(dat$G)),
                    logLambda = rep(log(lambda),length.out = 2),
                    thetaIn = matrix(0.0,nlevels(dat$G)-1,
                                     nlevels(proportionGroup)),
                    MIn = MIn,
                    tmixpIn = matrix(stats::qlogis(tMixture), n, nlevels(dat$G)),
                    logDf = matrix(log(tDf), n, nlevels(dat$G)),
                    kw1 = matrix(log(aParam), n, nlevels(dat$G)),
                    kw2 = matrix(log(bParam), n, nlevels(dat$G))
                    )
    }
    np <- Reduce("+",lapply(par,length))

    
##### Map for TMB #####
    ## variance parameter map
    if(independent){
        corparMap <- factor(rep(NA,length(par$corpar)))
    }else if(equalVariance){
        corparMap <- factor(row(par$corpar))
    }else{
        corparMap <- factor(rep(1:length(par$corpar),length = length(par$corpar)))
    }
    if(equalVariance){
        sigmaMap <- factor(row(par$logSigma))
        tMap <- factor(row(par$tmixpIn))
    }else{
        sigmaMap <- factor(1:length(par$logSigma))
        tMap <- factor(1:length(par$tmixpIn))
    }

    tMixMap <- tMap
    if(!estimateTMix)
        tMixMap[] <- NA
    tDfMap <- tMap
    if(!estimateTDf)
        tDfMap[] <- NA

    KWMapA <- tMap
    if(!estimateKWAParam)
        KWMapA[] <- NA

    ## logA and logB are not identifiable when the KW is scaled and moved
    KWMapB <- tMap
    if(!estimateKWBParam)
        KWMapB[] <- NA

    
    map <- list(logSigma = sigmaMap,
                corpar = corparMap,
                logLambda = factor(rep(1:length(lambda),length.out = length(par$logLambda))),
                MIn = factor(CMA_map),
                tmixpIn = factor(tMixMap),
                logDf = factor(tDfMap),
                kw1 = factor(KWMapA),
                kw2 = factor(KWMapB)
                )

    if(nrow(par$commonMu) > 0){
        cMuMap <- matrix(1:length(par$commonMu),nrow(par$commonMu),ncol(par$commonMu))
        cMuMap[1,] <- NA    
        map$commonMu <- factor(cMuMap)
    }
    if(all(is.na(map$commonMu))){
        map$logLambda <- factor(c(1,NA))
    }
    if(isTRUE(all.equal(lp_penalty,0))){
        map$logLambda <- factor(par$logLambda*NA)
        rnd <- c()
    }else{
        rnd <- c("mu","commonMu")
    }
    npfix <- Reduce("+",lapply(map,function(x)sum(is.na(x))))
    if(np-npfix > nrow(y) / 10)
        warning(sprintf("The model has %s parameters with only %s observations - %.2f observations per parameter. There is a high risk of overfitting.",np-npfix,nrow(y),nrow(y)/(np-npfix))) 
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
    opt$double_objective <- obj$env$f(obj$env$last.par.best, type="double")
    opt$gradient <- obj$gr(opt$par)
    
    rp <- obj$report(obj$env$last.par.best)

    muNames <- list(rownames(dat$covar),
                    colnames(y),
                    levels(group))
    commonMuNames <- list(rownames(dat$commonCovar),
                          colnames(y))

 
    if(doSdreport){
        opt$hessian <- stats::optimHess(opt$par, obj$fn, obj$gr)
        sdr <- TMB::sdreport(obj,opt$par, hessian.fixed = opt$hessian, getReportCovariance = getReportCovariance)
    }else{
        opt$hessian <- NULL
        sdr <- NULL
    }

    xlevels <- stats::.getXlevels(stats::terms(mf), mf)
    xlevelsCommon <- stats::.getXlevels(stats::terms(mfCommon), mfCommon)

    res <- list(call = cl,
                terms = stats::terms(mf),
                termsCommon = stats::terms(mfCommon),
                xlevels = xlevels,
                xlevelsCommon = xlevelsCommon,
                all_vars = colnames(stats::get_all_vars(stats::terms(mf), data)),
                all_varsCommon = colnames(stats::get_all_vars(stats::terms(mfCommon), data)),
                drop.unused.levels = drop.unused.levels,

                opt = opt,
                rp = rp,
                pl = obj$env$parList(opt$par),
                sdr = sdr,

                silent = silent,
                lp_penalty = lp_penalty,

                confusionMatrixArray = simplify2array(rp$Mvec),
                confusionLevelTypes = confusionLevelTypes,

                overfit = np-npfix > nrow(y) / 10,                
                muNames = muNames,
                commonMuNames = commonMuNames,

                y=y,
                group = group,
                data = data,

                proportionGroup = proportionGroup,
                confusionGroup = confusionGroup,
                
                tmb_data = obj$env$data,
                tmb_map = map,
                tmb_random = rnd,
                                
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
    
         
