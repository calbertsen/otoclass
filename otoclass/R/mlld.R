squeeze <- function(u, eps = .Machine$double.eps){
    (1.0 - eps) * (u - 0.5) + 0.5
}


##' Maximum likelihood linear discrimination
##'
##' @param proportionGroup 
##' @param confusionGroup 
##' @param ... 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats nlminb cor
##' @importFrom methods as
##' @importFrom TMB MakeADFun sdreport
##' @importFrom expm logm
##' @importFrom lme4 findbars nobars lFormula
##' @export
mlld <- function(## Data related
                 y,
                 group,
                 data = NULL,
                 formula = ~1,
                 formulaCommon = ~1,
                 formulaProportion = ~1,
                 formulaLogScale = ~ -1,
                 ##formulaDispersion = ~ -1,
                 ##proportionGroup = factor(rep("Baseline",nrow(y))),
                 ##confusionGroup = NULL,
                 ## dispersionGroup = factor(rep(NA,nrow(y))),
                 ## Penalty related
                 lp_penalty = NA,       # -1: Student's t penalty, 0: REML, 1: Lasso / Laplace prior, 2: Ridge / Gaussian prior, p: Lp penalty
                 lambda = 0.4,
                 estimateLambda = TRUE,
                 REML = FALSE,
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
                 doSdreport = TRUE,
                 getReportCovariance = FALSE,
                 equalVariance = TRUE,
                 confusionMatrixList = NULL,
                 groupConversionList = list(),
                 confusionLevelTypes = rep("Known",ifelse(is.null(ncol(group)),1,ncol(group))),  ## Possible values: c("Known","Fixed","Estimate")
                 observationType = c("MVMIX","SNP1","SNP2","SPL_AR1","FS1_AR1"),
                 forceMeanIncrease = FALSE,
                 ## na.pass = FALSE,
                 fixZeroGradient = TRUE,
                 SNP2dm = TRUE,
                 SPLknots = plogis(seq(qlogis(1/(ncol(y)+1)),qlogis(ncol(y)/(ncol(y)+1)),len=7)),
                 Nefd = 5,
                 lower = list(UlogSd = -10,
                              UThetalogSd = -10,
                              UComlogSd = -10,
                              betaTheta = -10),
                 upper = list(betaTheta = 10),
                 ...){
    
    cl <- match.call()
    observationType <- match.argInt(observationType) - 1

##### Checks #####
    if(observationType != 2){
        if(!is.matrix(y))
            if(is.vector(y)){
                y <- matrix(y,ncol = 1)
            }else{
                stop("y must be a matrix of observations")
            }
        yUse <- t(y)
    }else{
        if(!is.array(y))
            stop("y must be an array of allele observations with dimension #alleles x #loci x #individuals")
        yUse <- y
    }
    nObs <- ifelse(observationType==2, dim(yUse)[3L], ncol(yUse))
    ##nFeat <- ifelse(observationType==2, dim(yUse)[2L], nrow(yUse))
    nFeat <- switch(as.character(observationType),
                    "2" = dim(yUse)[2L],
                    "3" = length(SPLknots)+1,
                    "4" = Nefd * 2 + 1,
                    nrow(yUse))

    SNP2dm <- rep(SNP2dm,length.out = nFeat)
                       
    if(!is.factor(group) && !is.data.frame(group)){
        group <- factor(group)
    }

    if(is.data.frame(group)){
        gisf <- sapply(group, is.factor)
        if(any(!gisf)){
            for(i in which(gisf))
                group[[i]] <- factor(group[[i]])
        }
    }else{
        group <- data.frame(G = group)
    }

    groupLevels <- lapply(group, levels)
    groupNlevels <- sapply(group, nlevels)
    group <- do.call("cbind",lapply(group, function(x) as.integer(x) - 1))
    Ngroups <- unname(groupNlevels[1])
    
    if(length(groupConversionList) == 0){
        if(any(groupNlevels != Ngroups))
            stop("groupConversionList must be given when different number of groups are used")
        groupConversionList <- lapply(seq_along(groupNlevels), function(i) diag(1,Ngroups))
    }else if(length(groupConversionList) != ncol(group) - 1){
            stop("Length of groupConversionList must be one less than the number of columns of group.")
    }else{
        groupConversionList <- c(list(diag(1,Ngroups)), groupConversionList)
        if(!(all(sapply(groupConversionList,nrow) == groupNlevels) &&
             all(sapply(groupConversionList,ncol) == Ngroups))){
            stop("group conversion matrices have wrong dimensions.")
        }
    }
    
    ## if(!is.factor(proportionGroup))
    ##     proportionGroup <- factor(proportionGroup)
    ## if(!is.factor(confusionGroup))
    ##     confusionGroup <- factor(confusionGroup)
    ## if(!is.factor(dispersionGroup))
    ##     dispersionGroup <- factor(dispersionGroup)
    if(!(length(lambda) == 2 || length(lambda)==1) | any(lambda <= 0))
        stop("Lambda must be a positive scalar.")
    ## if((!identical(formula, ~1) | !identical(formulaCommon, ~1)) & is.null(data))
    ##     stop("When a formula is specified, data must be given.")
    if(!is.null(data) & !identical(nObs, nrow(data)))
        stop("data and y must have the same number of rows.")

    confusionLevelTypes <- c("Known","Fixed","Estimate")[pmatch(confusionLevelTypes,c("Known","Fixed","Estimate"), duplicates.ok = TRUE)]
    if(any(is.na(confusionLevelTypes)))
        stop('confusionLevelTypes must be a vector of "Known","Fixed","Estimate". Partial matches are allowed.')
    if((length(confusionLevelTypes) != ncol(group)))
        stop("The length of confusionLevelTypes must match the number of columns of group.")
    if(!is.null(confusionMatrixList) && !((length(confusionMatrixList) == sum(confusionLevelTypes=="Fixed")) || (length(confusionMatrixList) == length(confusionLevelTypes))))
        stop("...")

##### Prepare confusion matrix array #####
    hasAll <- !is.null(confusionMatrixList) && (length(confusionMatrixList) == length(confusionLevelTypes))
    CMA <- lapply(groupNlevels, function(ii) diag(1,ii))
    
    if(hasAll){
        CMA <- confusionMatrixList
    }else{
        for(i in seq_along(confusionLevelTypes)){
            if(confusionLevelTypes[i] == "Fixed"){
                isFixed <- confusionLevelTypes[1:i]=="Fixed"
                CMA[[i]] <- confusionMatrixList[[sum(isFixed)]]
            }else if(confusionLevelTypes[i] == "Estimate"){
                CMA[[i]][] <- 0.1/nrow(CMA[[i]])
                diag(CMA[[i]]) <- 0.9
            }
        }
    }

    CMA_mapL <- relist(seq_along(unlist(CMA)), CMA)
    for(i in seq_along(CMA_mapL)){
        if(confusionLevelTypes[i] %in% c("Estimate")){
            diag(CMA_mapL[[i]]) <- NA
        }else{
            CMA_mapL[[i]][] <- NA
        }
    }
    CMA_map <- factor(unlist(CMA_mapL))
    
  
##### Prepare model matrix #####
    lc <- lme4::lmerControl(check.nobs.vs.rankZ = "ignore",
                            check.nobs.vs.nlev = "ignore",
                            check.nlev.gtreq.5 = "ignore",
                            check.nlev.gtr.1 = "ignore",
                            check.nobs.vs.nRE= "ignore",
                            check.rankX = "ignore",
         check.scaleX = "ignore",
         check.formula.LHS = "ignore")
    
    if(is.null(data) & identical(formula,~1) & identical(formulaCommon, ~1))
        data <- data.frame(ID = 1:nObs)
    if(!is.null(data) & !is.data.frame(data))
        data <- as.data.frame(data)

    mf <- model.frame(lme4::nobars(formula), data,
                              na.action = na.pass)

    X <- Matrix::sparse.model.matrix(terms(mf),
                                     data = mf,
                                     transpose = TRUE,
                                     row.names = FALSE,
                                     drop.unused.levels = drop.unused.levels)
    if(inherits(X,"dgCMatrix"))
        X <- as(X,"dgTMatrix")    

    if(is.null(lme4::findbars(formula))){
        Z <- list()
        U <- array(0,dim = c(0))
        attr(U,"rdim") <- integer(0)
        attr(U,"cdim") <- integer(0)
        attr(U,"adim") <- integer(0)
        Ucor <- array(0, dim = c(0))
        attr(Ucor,"rdim") <- integer(0)
        attr(Ucor,"cdim") <- integer(0)
        UlogSd <- array(0, dim = c(0))
        attr(UlogSd,"rdim") <- integer(0)
        attr(UlogSd,"cdim") <- integer(0)        
    }else{
        rtZ <- lme4::lFormula(formula,data, na.action = na.pass, control = lc)$reTrms
        Z <- lapply(rtZ$Ztlist,function(xx){
            as(xx,"dgTMatrix")
        })       
        Znms <- rtZ$cnms
        Zrdim <- sapply(rtZ$cnms,length)
        Zcdim <- sapply(rtZ$flist,nlevels)
        U <- array(numeric(sum(Zrdim * Zcdim * nFeat * Ngroups)))
        attr(U,"rdim") <- as.integer(Zrdim * Zcdim)
        attr(U,"cdim") <- as.integer(rep(nFeat, length(Z)))
        attr(U,"adim") <- as.integer(rep(Ngroups, length(Z)))        
        n <- Zrdim * (Zrdim - 1) / 2
        Ucor <- array(0, dim = sum(n * nFeat * Ngroups))
        attr(Ucor,"rdim") <- as.integer(n)
        attr(Ucor,"cdim") <- as.integer(rep(nFeat, length(Z)))
        attr(Ucor,"adim") <- as.integer(rep(Ngroups, length(Z)))
        UlogSd <- array(2, dim = sum(Zrdim * nFeat * Ngroups))
        attr(UlogSd,"rdim") <- as.integer(Zrdim)
        attr(UlogSd,"cdim") <- as.integer(rep(nFeat, length(Z)))
        attr(UlogSd,"adim") <- as.integer(rep(Ngroups, length(Z)))
    }

##### Prepare common model matrix #####
    formulaCommon <- update.formula(formulaCommon, ~ . +1)
    mfCommon <- model.frame(lme4::nobars(formulaCommon), data,
                              na.action = na.pass)
    XCom <- Matrix::sparse.model.matrix(terms(mfCommon),
                                        data = mfCommon,
                                        transpose = TRUE,
                                        row.names = FALSE,
                                        drop.unused.levels = drop.unused.levels)
    if(inherits(XCom,"dgCMatrix"))
        XCom <- as(XCom,"dgTMatrix")

    if(is.null(lme4::findbars(formulaCommon))){
        ZCom <- list()
        UC <- array(0,dim = c(0))
        attr(UC,"rdim") <- integer(0)
        attr(UC,"cdim") <- integer(0)
        UCcor <- array(0, dim = c(0))
        attr(UCcor,"rdim") <- integer(0)
        attr(UCcor,"cdim") <- integer(0)
        UClogSd <- array(0, dim = c(0))
        attr(UClogSd,"rdim") <- integer(0)
        attr(UClogSd,"cdim") <- integer(0)
   }else{
        rtZC <- lme4::lFormula(formulaCommon,data, na.action = na.pass, control = lc)$reTrms
        ZCom <- lapply(rtZC$Ztlist,function(xx){
            as(xx,"dgTMatrix")
        })
        ZCnms <- rtZC$cnms
        ZCrdim <- sapply(rtZC$cnms,length)
        ZCcdim <- sapply(rtZC$flist,nlevels)
        UC <- array(numeric(sum(ZCrdim * ZCcdim * nFeat)))
        attr(UC,"rdim") <- as.integer(ZCrdim * ZCcdim)
        attr(UC,"cdim") <- as.integer(rep(nFeat, length(ZCom)))
        n <- ZCrdim * (ZCrdim - 1) / 2
        UCcor <- array(0, dim = sum(n*nFeat))
        attr(UCcor,"rdim") <- as.integer(n)
        attr(UCcor,"cdim") <- as.integer(rep(nFeat,length(ZCom)))
        UClogSd <- array(2, dim = sum(ZCrdim*nFeat))
        attr(UClogSd,"rdim") <- as.integer(ZCrdim)
        attr(UClogSd,"cdim") <- as.integer(rep(nFeat,length(ZCom)))
    }

    
##### Prepare proportion model matrix #####
    formulaLogScale <- update.formula(formulaLogScale, ~ . +1)
    mfLogScale <- model.frame(lme4::nobars(formulaLogScale), data,
                              na.action = na.pass)
    XLogScale <- Matrix::sparse.model.matrix(terms(mfLogScale),
                                        data = mfLogScale,
                                        transpose = TRUE,
                                        row.names = FALSE,
                                        drop.unused.levels = drop.unused.levels)
    if(inherits(XLogScale,"dgCMatrix"))
        XLogScale <- as(XLogScale,"dgTMatrix")

    if(!is.null(lme4::findbars(formulaLogScale)))
        warning("formulaLogScale does not allow random effets. Random effect terms were removed.")
    
##### Prepare proportion model matrix #####
    mfTheta <- model.frame(lme4::nobars(formulaProportion), data,
                              na.action = na.pass)
    XTheta <- Matrix::sparse.model.matrix(terms(mfTheta),
                                        data = mfTheta,
                                        transpose = TRUE,
                                        row.names = FALSE,
                                        drop.unused.levels = drop.unused.levels)
    if(inherits(XTheta,"dgCMatrix"))
        XTheta <- as(XTheta,"dgTMatrix")

    if(is.null(lme4::findbars(formulaProportion))){
        ZT <- list()
        UT <- array(0,dim = c(0))
        attr(UT,"rdim") <- integer(0)
        attr(UT,"cdim") <- integer(0)
        attr(UT,"adim") <- integer(0)
        UTcor <- array(0, dim = c(0))
        attr(UTcor,"rdim") <- integer(0)
        attr(UTcor,"cdim") <- integer(0)        
        UTlogSd <- array(0, dim = c(0))
        attr(UTlogSd,"rdim") <- integer(0)
        attr(UTlogSd,"cdim") <- integer(0)        
    }else{
        rtZT <- lme4::lFormula(formulaProportion,data, na.action = na.pass, control = lc)$reTrms
        ZT <- lapply(rtZT$Ztlist,function(xx){
            as(xx,"dgTMatrix")
        })
        ZTnms <- rtZT$cnms
        ZTrdim <- sapply(rtZT$cnms,length)
        ZTcdim <- sapply(rtZT$flist,nlevels)
        UT <- array(numeric(sum(ZTrdim * ZTcdim * (Ngroups-1))))
        attr(UT,"rdim") <- as.integer(ZTrdim * ZTcdim)
        attr(UT,"cdim") <- as.integer(rep(Ngroups-1, length(ZT)))
        n <- ZTrdim * (ZTrdim - 1) / 2
        UTcor <- array(0, dim = sum(n*(Ngroups-1)))
        attr(UTcor,"rdim") <- as.integer(n)
        attr(UTcor,"cdim") <- as.integer(rep(Ngroups-1, length(ZT)))
        UTlogSd <- array(2, dim = sum(ZTrdim * (Ngroups-1)))
        attr(UTlogSd,"rdim") <- as.integer(ZTrdim)
        attr(UTlogSd,"cdim") <- as.integer(rep(Ngroups-1, length(ZT)))
   }

    
##### Prepare proportion overdispersion model matrix #####
    ## mfDisp <- model.frame(formulaDispersion, data,
    ##                           na.action = na.fail)
    ## XDisp <- Matrix::sparse.model.matrix(terms(mfDisp),
    ##                                       data = mfDisp,
    ##                                       transpose = TRUE,
    ##                                       row.names = FALSE,
    ##                                       drop.unused.levels = drop.unused.levels)
    ## if(inherits(XDisp,"dgCMatrix"))
    ##     XDisp <- as(XDisp,"dgTMatrix")

##### Prepare confusion matrices #####
    MIn <- CMA
    for(i in seq_along(MIn)){
        if(dim(MIn[[1]])[1] == 1){
            MIn[[i]] <- 1
        }else{
            MIn[[i]] <- expm::logm(MIn[[i]] + diag(1e-10, dim(MIn[[i]])[1]))
            diag(MIn[[i]]) <- 1
        }
        MIn[[i]] <- log(MIn[[i]])
    }

    MInCMOE <- as.array(unname(unlist(MIn)))
    attr(MInCMOE,"rdim") <- as.integer(unname(sapply(MIn,nrow)))
    attr(MInCMOE,"cdim") <- as.integer(unname(sapply(MIn,ncol)))

##### Handle NA in group #####
    ## Guse <- group
    ## if(any(is.na(Guse))){
    ##     indx <- which(is.na(Guse))
    ##     if(any(confusionLevelTypes[confusionGroup[indx]] != "Unknown"))
    ##         stop("All observations with NA in group must have a confusionGroup with confusionLevelType equal to 'Unknown'")
    ##     Guse[indx] <- levels(group)[1]
    ## }

    ## if(any(is.na(dispersionGroup)))
    ##     stop("dispersionGroup can not have NA")
    
    ## if(!is.factor(dispersionGroup))
    ##     dispersionGroup <- factor(dispersionGroup)

#####

    increaseFirstCoordinate <- ifelse(forceMeanIncrease || all(is.na(group)), 1L, 0L)

    xx <- lapply(groupConversionList, function(x) t(x)%*%x)
    xxIndx <- apply(group,2,function(x) !all(is.na(x)))
    if(sum(xxIndx) > 0){
        xxR <- Reduce("+", xx[xxIndx])
        ee <- eigen(xxR)
    }else{
        ee <- list(value = 1)
    }
    if(any(abs(ee$value) < 1e-8)){
        eeV <- ee$vector
        identifyMatrix <- unique(do.call("rbind",lapply(as.data.frame(eeV[,ee$value < 1e-8,drop=FALSE]),
                                                        function(x){ i <- which(abs(x) > 1e-4);
                                                            as.matrix(subset(expand.grid(i,i),Var1 < Var2)) - 1})
                                         ))
        identifyMatrix <- identifyMatrix[order(identifyMatrix[,1]),,drop=FALSE]

    }else{
        identifyMatrix <- matrix(NA_integer_,0,2)
    }

    
##### Data for TMB #####
    dat <- list(model = observationType,
                G = group,
                Gnlevels = groupNlevels,
                Gconversion = groupConversionList,
                Y = yUse,
                X = X,
                XCom = XCom,
                XLogScale = XLogScale,
                XTheta = XTheta,
                Z = Z,
                ZCom = ZCom,
                ZTheta = ZT,
                ## proportionGroup = proportionGroup,
                ## confusionGroup = confusionGroup,
                ## dispersionGroup = dispersionGroup,
                penalty = ifelse(is.na(lp_penalty),0,lp_penalty),
                Y_pred = array(0,dim=c(0,0)),
                G_pred = matrix(0L, 0, 0),
                X_pred = as(matrix(0,0,0),"dgTMatrix"),
                XCom_pred = as(matrix(0,0,0),"dgTMatrix"),
                XLogScale_pred = as(matrix(0,0,0),"dgTMatrix"),
                XTheta_pred = as(matrix(0,0,0),"dgTMatrix"),
                Z_pred = list(),
                ZCom_pred = list(),
                ZTheta_pred = list(),
                ## proportionGroup_pred = factor(),
                ## dispersionGroup_pred = factor(),
                increaseFirstCoordinate = increaseFirstCoordinate,
                identifyMatrix = identifyMatrix,
                knots = SPLknots
                )

##### Parameters for TMB #####
    n <- nFeat ## nrow(dat$Y)
    ## par list
    mndim <- c(nrow(dat$X), nFeat * ifelse(observationType==2,nrow(dat$Y)-1,1), dat$Gnlevels[1])
    mn <- array(rnorm(prod(mndim),0,0.001), dim = mndim) ## 
    if(observationType == 2){ ## Better starting values for SNP2
        sv <- sapply(seq_len(dim(dat$Y)[2]), function(i)
            unlist(lapply(split(as.data.frame(t(dat$Y[,i,])),factor(dat$G[,1]+1,seq_len(dat$Gnlevels[1]))),function(x){ v <- colSums(as.matrix(x)); qlogis(squeeze(v[1]/sum(v))) }))
            )        
        mn[1,,][which(!is.nan(sv))] <- unname(sv[which(!is.nan(sv))])
        if(forceMeanIncrease){
            isv <- !is.nan(sv[,1])
            sv[,1] <- approx(x = c(0,seq_len(dat$Gnlevels[1])[isv],dat$Gnlevels[1]+1),
                             y = c(-10,sv[isv,1],10), xout = seq_len(dat$Gnlevels[1]))$y
            if(!all(diff(sv[,1]) > 0))
                warning("Mean values for the first feature does not apear to increase with groups.")
            mn[1,1,] <- unname(c(sv[1,1], log(diff(sv[,1]))))
        }
    }
    ##mn <- sapply(levels(dat$G),function(i)apply(dat$X[,dat$G==i,drop=FALSE],1,mean))
    if(independent || observationType %in% c(1,2,3,4)){
        corcalc <- 0
    }else{
        ##corcalc <- t(chol(stats::cor(t(dat$Y))))[lower.tri(t(chol(stats::cor(t(dat$Y)))),diag = FALSE)]
        corcalc <- 0
    }

    if(observationType %in% c(1,2)){
        logSd <- numeric(0)
    }else if(observationType %in% c(3,4)){
        logSd <- numeric(2)
    }else{
        logSd <- log(apply(dat$Y,1,sd))+2
    }
    
    
    if(!is.null(cl$parlist)){
        par <- eval(cl$parlist)
        if(!is.null(cl$lambda))
            par$logLambda <- log(lambda)
    }else{
        par <- list(mu = mn,
                    commonMu = matrix(0,nrow(dat$XCom),nFeat * ifelse(observationType==2,nrow(dat$Y)-1,1)),
                    logSigma = matrix(logSd,
                                      ifelse(observationType%in%c(1,2),0,ifelse(observationType%in%c(3,4),2,n)),
                                      dat$Gnlevels[1]),
                    corpar = matrix(corcalc,
                                    ifelse(observationType%in%c(1,2,3,4),0,(n*n-n)/2),
                                    dat$Gnlevels[1]),
                    U = U,
                    Ucor = Ucor,
                    UlogSd = UlogSd,
                    UCom = UC,
                    UComcor = UCcor,
                    UComlogSd = UClogSd,
                    UTheta = UT,
                    UThetacor = UTcor,
                    UThetalogSd = UTlogSd,
                    logLambda = rep(log(lambda),length.out = 2),
                    betaLogScale = array(0, dim = c(nrow(dat$XLogScale), nFeat, Ngroups)),
                    betaTheta = matrix(0, nrow(dat$XTheta), Ngroups-1),
                    ## thetaIn = matrix(0.0,groupNlevels[1]-1,
                    ##                  nlevels(proportionGroup)),
                    MIn = MInCMOE,
                    tmixpIn = matrix(qlogis(tMixture+0.01*as.numeric(estimateTMix)), n, dat$Gnlevels[1]),
                    logDf = matrix(log(tDf), n, dat$Gnlevels[1])
                    )
    }
    
    np <- Reduce("+",lapply(par[!(names(par) %in% c("U","UCom","UTheta"))],length))

    
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

    ## if(nrow(dat$XDisp) == 0){
    ##     map$logSdDispersion = factor(NA)
    ##     map$dispersion = factor(rep(NA, length(par$dispersion)))        
    ## }

    if(nrow(par$commonMu) > 0){
        cMuMap <- matrix(1:length(par$commonMu),nrow(par$commonMu),ncol(par$commonMu))
        cMuMap[1,] <- NA    
        map$commonMu <- factor(cMuMap)
    }
    if(nrow(par$betaLogScale) > 0){
            blsMap <- array(1:length(par$betaLogScale),dim = dim(par$betaLogScale))
        if(observationType %in% c(1,3,4)){
            blsMap[] <- NA
        }else if(observationType == 0){
            blsMap[1,,] <- NA
        }else if(observationType == 2){ # For observationType == 2, it is used for dirichlet multinomial
            blsMap[,!SNP2dm,] <- NA
            par$betaLogScale[,!SNP2dm,] <- Inf
        }                        
        map$betaLogScale <- factor(blsMap)
    }
    if(all(is.na(map$commonMu))){
        map$logLambda <- factor(c(1,NA))
    }
    if(!REML){
        if(all(is.na(lp_penalty))){
            map$logLambda <- factor(par$logLambda*NA)
            rnd <- c()
        }else{
            rnd <- c("mu","commonMu")
        }
        RE <- c("U","UCom","UTheta")
    }else{
        if(all(is.na(lp_penalty)))
            map$logLambda <- factor(par$logLambda*NA)
        rnd <- c()
        RE <- c("mu","commonMu","U","UCom","UTheta")

    }
    npfix <- Reduce("+",lapply(map,function(x)sum(is.na(x))))
    if(np-npfix > nObs / 10)
        warning(sprintf("The model has %s parameters with only %s observations - %.2f observations per parameter. There is a high risk of overfitting.",np-npfix,nObs,nObs/(np-npfix)))
#### return(list(dat=dat,par=par,map=map))
##### Make TMB object #####
    obj <- TMB::MakeADFun(dat,par,map,
                          silent = silent, profile = rnd,
                          random = RE,
                          DLL = "otoclass")

    if(fixZeroGradient && any((ogr <- obj$gr()) == 0)){
        toFix <- split(ogr == 0,names(obj$par))
        for(nn in names(toFix))
            if(any(toFix[[nn]])){
                if(is.null(map[[nn]]))
                    map[[nn]] <- factor(seq_along(toFix[[nn]]))
                map[[nn]][!is.na(map[[nn]])][toFix[[nn]]] <- NA
                map[[nn]] <- factor(map[[nn]])
            }
        obj <- TMB::MakeADFun(dat,par,map,
                              silent = silent, profile = rnd,
                              random = RE,
                              DLL = "otoclass")
    }

    
    if(onlyObj)
        return(obj)

    low <- rep(-Inf, length(obj$par))
    for(nn in names(lower))
        low[names(obj$par) %in% nn] <- as.vector(lower[[nn]])
    upp <- rep(Inf, length(obj$par))
    for(nn in names(upper))
        upp[names(obj$par) %in% nn] <- as.vector(upper[[nn]])
    
    opt <- stats::nlminb(obj$par,
                         obj$fn,
                         obj$gr,
                         lower = low,
                         upper = upp,
                         control = control)
    opt$double_objective <- obj$env$f(obj$env$last.par.best, type="double")
    rp <- obj$report(obj$env$last.par.best)

    muNames <- list(rownames(dat$X),
                    colnames(y),
                    groupLevels[[1]])
    commonMuNames <- list(rownames(dat$XCom),
                          colnames(y))

 
    if(doSdreport){
        sdr <- TMB::sdreport(obj,opt$par, getReportCovariance = getReportCovariance)
    }else{
        sdr <- NULL
    }

    xlevels <- .getXlevels(terms(mf), mf)
    xlevelsCommon <- .getXlevels(terms(mfCommon), mfCommon)
    xlevelsTheta <- .getXlevels(terms(mfTheta), mfTheta)
    ## xlevelsDisp <- .getXlevels(terms(mfDisp), mfDisp)

    res <- list(call = cl,
                terms = terms(mf),
                termsCommon = terms(mfCommon),
                termsTheta = terms(mfTheta),
                ## termsDisp = terms(mfDisp),
                xlevels = xlevels,
                xlevelsCommon = xlevelsCommon,
                xlevelsTheta = xlevelsTheta,
                ## xlevelsDisp = xlevelsDisp,
                all_vars = colnames(get_all_vars(terms(mf), data)),
                all_varsCommon = colnames(get_all_vars(terms(mfCommon), data)),
                all_varsTheta = colnames(get_all_vars(terms(mfTheta), data)),
                ## all_varsDisp = colnames(get_all_vars(terms(mfDisp), data)),
                drop.unused.levels = drop.unused.levels,

                opt = opt,
                rp = rp,
                pl = obj$env$parList(opt$par),
                sdr = sdr,

                silent = silent,
                lp_penalty = lp_penalty,

                confusionMatrixList = (rp$Mvec),
                confusionLevelTypes = confusionLevelTypes,

                overfit = np-npfix > nObs / 10,                
                muNames = muNames,
                commonMuNames = commonMuNames,

                y=y,
                group = group,
                groupLevels = groupLevels,
                
                data = data,

                ## proportionGroup = proportionGroup,
                ## confusionGroup = confusionGroup,
                
                tmb_data = obj$env$data,
                tmb_map = map,
                tmb_random = rnd,

                identifyMatrix = identifyMatrix,
                                
                varBetweenGroups = NA,
                varWithinGroups = NA
                )
                
    class(res) <- "mlld"
    if(observationType == 0){
        res$varBetweenGroups <- calculateVarBetweenGroups(res)
        res$varWithinGroups <- calculateVarWithinGroups(res)
    }
    
    return(res)    
}

##' @export
print.mlld <- function(x, ...){
    print(x$opt)
}
    



##' @export
AIC.mlld <- function(x, k = 2, ...){
    k * x$opt$objective + k * length(x$opt$par)
}
