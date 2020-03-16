##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param proportionGroup 
##' @param confusionGroup 
##' @param ... 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats nlminb cor
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
                 dispersionGroup = factor(rep(1,nrow(y))),
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
                CMA[,,i] <- diag(1,nrow = dim(CMA)[1])
            }else if(confusionLevelTypes[i] == "Fixed"){
                isFixed <- confusionLevelTypes[1:i]=="Fixed"
                CMA[,,i] <- confusionMatrixArray[,,sum(isFixed)]
            }else if(confusionLevelTypes[i] == "Estimate"){
                tmpMat <- matrix(0.1/dim(CMA)[1],nrow = dim(CMA)[1], ncol = dim(CMA)[2])
                diag(tmpMat) <- 0.9
                CMA[,,i] <- tmpMat
            }else if(confusionLevelTypes[i] == "Unknown"){
                CMA[,,i] <- 1/dim(CMA)[1]
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

    mf <- model.frame(formula, data)

    X <- Matrix::sparse.model.matrix(terms(mf),
                                     data = mf,
                                     transpose = TRUE,
                                     row.names = FALSE,
                                     drop.unused.levels = drop.unused.levels)
    if(inherits(X,"dgCMatrix"))
        X <- as(X,"dgTMatrix")

##### Prepare common model matrix #####
    formulaCommon <- update.formula(formulaCommon, ~ . +1)
    mfCommon <- model.frame(formulaCommon, data)
    XCom <- Matrix::sparse.model.matrix(terms(mfCommon),
                                        data = mfCommon,
                                        transpose = TRUE,
                                        row.names = FALSE,
                                        drop.unused.levels = drop.unused.levels)
    if(inherits(XCom,"dgCMatrix"))
        XCom <- as(XCom,"dgTMatrix")

##### Prepare confusion matrices #####
    MIn <- CMA
    for(i in 1:dim(MIn)[3]){
        if(dim(MIn)[1] == 1){
            MIn[,,i] <- 1
        }else{
            MIn[,,i] <- expm::logm(MIn[,,i])
            diag(MIn[,,i]) <- 1
        }
        MIn[,,i] <- log(MIn[,,i])
    }

##### Handle NA in group #####
    Guse <- group
    if(any(is.na(Guse))){
        indx <- which(is.na(Guse))
        if(any(confusionLevelTypes[confusionGroup[indx]] != "Unknown"))
            stop("All observations with NA in group must have a confusionGroup with confusionLevelType equal to 'Unknown'")
        Guse[indx] <- levels(group)[1]
    }

    if(any(is.na(dispersionGroup)))
        stop("dispersionGroup can not have NA")
    
    if(!is.factor(dispersionGroup))
        dispersionGroup <- factor(dispersionGroup)
    
##### Data for TMB #####
    dat <- list(model = ifelse(independent,0L,1L),
                G = Guse,
                Y = t(y),
                X = X,
                XCom = XCom,
                proportionGroup = proportionGroup,
                confusionGroup = confusionGroup,
                dispersionGroup = dispersionGroup,
                penalty = lp_penalty,
                Y_pred = array(0,dim=c(0,0)),
                X_pred = as(matrix(0,0,0),"dgTMatrix"),
                XCom_pred = as(matrix(0,0,0),"dgTMatrix"),
                proportionGroup_pred = factor(),
                dispersionGroup_pred = factor()
                )
    
##### Parameters for TMB #####
    n <- nrow(dat$Y)
    ## par list
    mndim <- c(nrow(dat$X), nrow(dat$Y), nlevels(dat$G))
    mn <- array(rnorm(prod(mndim),0,0.001), dim = mndim)
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
        par <- list(mu = mn,
                    commonMu = matrix(0,nrow(dat$XCom),nrow(dat$Y)),
                    logSigma = matrix(log(apply(dat$Y,1,sd))+2,n,nlevels(dat$G)),
                    corpar = matrix(corcalc,(n*n-n)/2,nlevels(dat$G)),
                    logLambda = rep(log(lambda),length.out = 2),
                    thetaIn = matrix(0.0,nlevels(dat$G)-1,
                                     nlevels(proportionGroup)),
                    MIn = MIn,
                    tmixpIn = matrix(qlogis(tMixture+0.01*as.numeric(estimateTMix)), n, nlevels(dat$G)),
                    logDf = matrix(log(tDf), n, nlevels(dat$G)),
                    logSdDispersion = 0,
                    dispersion = matrix(0.0,nlevels(dat$G)-1,
                                        nlevels(dispersionGroup))
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

    
    
    map <- list(logSigma = sigmaMap,
                corpar = corparMap,
                logLambda = factor(rep(1:length(lambda),length.out = length(par$logLambda))),
                MIn = factor(CMA_map),
                tmixpIn = factor(tMixMap),
                logDf = factor(tDfMap)
                )

    if(nlevels(dispersionGroup) == 1){
        map$logSdDispersion = factor(NA)
        map$dispersion = factor(rep(NA, length(par$dispersion)))        
    }

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
    if(nlevels(dispersionGroup) > 1){
        rnd <- c(rnd,"dispersion")
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
    low <- rep(-Inf, length(obj$par))
    low[names(obj$par) %in% "logSdDispersion"] <- -10
    opt <- stats::nlminb(obj$par,
                         obj$fn,
                         obj$gr,
                         lower = low,
                         control = control)
    opt$double_objective <- obj$env$f(obj$env$last.par.best, type="double")
    rp <- obj$report(obj$env$last.par.best)

    muNames <- list(rownames(dat$X),
                    colnames(y),
                    levels(group))
    commonMuNames <- list(rownames(dat$XCom),
                          colnames(y))

 
    if(doSdreport){
        sdr <- TMB::sdreport(obj,opt$par, getReportCovariance = getReportCovariance)
    }else{
        sdr <- NULL
    }

    xlevels <- .getXlevels(terms(mf), mf)
    xlevelsCommon <- .getXlevels(terms(mfCommon), mfCommon)

    res <- list(call = cl,
                terms = terms(mf),
                termsCommon = terms(mfCommon),
                xlevels = xlevels,
                xlevelsCommon = xlevelsCommon,
                all_vars = colnames(get_all_vars(terms(mf), data)),
                all_varsCommon = colnames(get_all_vars(terms(mfCommon), data)),
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
    
         
