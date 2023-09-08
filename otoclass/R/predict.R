getProportion <- function(f, data, indx){
    if(missing(indx))
        indx <- seq_along(f$muNames[[3]])
    if(missing(data))
        data <- f$data
    if(is.null(f$sdr)){
        betaTheta <- as.vector(f$pl$betaTheta)
        Sigm <- matrix(NA_real_,length(betaTheta),length(betaTheta))
    }else{
        p0 <- as.list(fit$sdr,"Estimate")$betaTheta
        if(!is.null(fit$tmb_map$betaTheta)){
            m0 <- fit$tmb_map$betaTheta
        }else{
            m0 <- seq_along(p0)
        }
        Sigm <- matrix(0,length(p0),length(p0))
        ssdr <- summary(f$sdr)
        betaTheta <- p0[] # unname(ssdr[rownames(ssdr)=="betaTheta",1])
        Sigm[!is.na(m0),!is.na(m0)] <- f$sdr$cov.fixed[rownames(f$sdr$cov.fixed) == "betaTheta",colnames(f$sdr$cov.fixed) == "betaTheta"]
    }
    X <- model.matrix(f$termsTheta,data = data, xlev=f$xlevelsTheta)
    fn <- function(p) {
        P <- matrix(p,nrow(f$pl$betaTheta))
        x <- cbind(X %*% P,0)
        qlogis(apply(x, 1, function(y) exp(y) / sum(exp(y)))[indx,])
    }
    rEst <- unname(fn(betaTheta))
    g <- numDeriv::jacobian(fn, betaTheta)
    rVar <- (g %*% Sigm %*% t(g))
    rSd <- sqrt(diag(rVar))
    dd <- list(Estimate = plogis(rEst),
               CIlow = plogis(rEst - 2 * rSd),
               CIhigh = plogis(rEst + 2 * rSd))
    dd <- lapply(dd, function(xx){
        if(is.matrix(xx))
            return(xx)
        return(matrix(xx, length(f$muNames[[3]]), nrow(data)))
    })
    attr(dd,"rEst") <- rEst
    attr(dd,"rVar") <- rVar
    rownames(dd[[1]]) <- rownames(dd[[2]]) <- rownames(dd[[3]]) <- f$muNames[[3]]
    colnames(dd[[1]]) <- colnames(dd[[2]]) <- colnames(dd[[3]]) <- rownames(data)
    dd <- lapply(dd, function(xx) xx[indx,,drop=FALSE])
    if(length(indx) == 1 || nrow(data) == 1){
        dnm <- paste(rownames(dd[[1]])[row(dd[[1]])], colnames(dd[[1]])[col(dd[[1]])], sep = "_")
        dd <- as.data.frame(lapply(dd,as.vector))
        rownames(dd) <- dnm
    }
    dd
}

##' Calculate group proportions for a data frame of covariates
##'
##' @param f 
##' @param data 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @export
getGroupProportion <- function(f, data,...){
    UseMethod("getGroupProportion")
}

##' @rdname getGroupProportion
##' @method getGroupProportion mlld
##' @export
getGroupProportion.mlld <- function(f, data, randEff = TRUE,...){
    require(lme4)
   if(missing(data))
        data <- f$data
    tmb_data <- f$tmb_data
    tmb_par <- f$pl
    tmb_map <- f$tmb_map
    tmb_random <- f$tmb_random
    tmb_profile <- f$tmb_profile

    ## Handle fixed effects
    
    tmb_data$XTheta_pred <- Matrix::sparse.model.matrix(f$termsTheta,
                                                   data = data,
                                                   transpose = TRUE,
                                                   row.names = FALSE,
                                                   xlev = f$xlevelsTheta)
    tmb_data$XTheta_pred <- as(tmb_data$XTheta_pred,"TsparseMatrix")
     ## Handle random effects
    if(!is.null(lme4::findbars(f$call$formulaProportion)) & randEff){
        mf <- model.frame(lme4::subbars(f$call$formulaProportion),f$data)
        xlvs <- .getXlevels(terms(mf), mf)
        dataSafe <- data
        xNA <- lapply(dataSafe[names(xlvs)], is.na)
        for(nn in names(xNA))
            dataSafe[[nn]][xNA[[nn]]] <- xlvs[[nn]][1]
        mf2 <- model.frame(terms(mf),data=dataSafe,xlev=.getXlevels(terms(mf), mf),na.action=na.pass)
        rtZT <- mkReTrms(findbars(f$call$formulaProportion),mf2,drop.unused.levels=FALSE)
        ZTtmp <- rtZT$Ztlist
        grpnames <- strsplit(gsub(".+\\| ","",names(ZTtmp)), ":")
        for(ii in seq_along(grpnames)){
            jj <- na.omit(match(grpnames[[ii]],names(xNA)))
            if(length(jj) == 0)
                next;
            kk <- apply(do.call("cbind",xNA[jj]),1,any)
            if(sum(kk) > 0)
                ZTtmp[[ii]][,kk] <- 0
        }    
        ZT <- lapply(ZTtmp,function(xx){
            as(xx,"TsparseMatrix")
        })
        ## ZTnms <- rtZT$cnms
        ## ZTrdim <- sapply(rtZT$cnms,length)
        ## ZTcdim <- sapply(rtZT$flist,nlevels)
        ## UT <- array(numeric(sum(ZTrdim * ZTcdim * (Ngroups-1))))
        ## attr(UT,"rdim") <- as.integer(ZTrdim * ZTcdim)
        ## attr(UT,"cdim") <- as.integer(rep(Ngroups-1, length(ZT)))
        ## n <- ZTrdim * (ZTrdim - 1) / 2
        ## UTcor <- array(0, dim = sum(n*(Ngroups-1)))
        ## attr(UTcor,"rdim") <- as.integer(n)
        ## attr(UTcor,"cdim") <- as.integer(rep(Ngroups-1, length(ZT)))
        ## UTlogSd <- array(2, dim = sum(ZTrdim * (Ngroups-1)))
        ## attr(UTlogSd,"rdim") <- as.integer(ZTrdim)
        ## attr(UTlogSd,"cdim") <- as.integer(rep(Ngroups-1, length(ZT)))
        tmb_data$ZTheta_pred <- ZT
   }

    
    obj <- TMB::MakeADFun(data = tmb_data,
                     parameters = tmb_par,
                     map = tmb_map,
                     random = tmb_random,
                     profile = tmb_profile,
                     silent = TRUE,
                     DLL = "otoclass")
    rp <- obj$report()
    sdr <- TMB::sdreport(obj, obj$par, f$opt$he)
    ssdr <- TMB::summary.sdreport(sdr)
    ##
    grpP <- ssdr[rownames(ssdr) == "GROUP_logitprobability",]
    fit_P <- plogis(do.call("structure",c(list(.Data=unname(grpP[,1])), attributes(rp$GROUP_logitprobability))))
    CIlow_P <- plogis(do.call("structure",c(list(.Data=unname(grpP[,1] - 2 * grpP[,2])), attributes(rp$GROUP_logitprobability))))
    CIhigh_P <- plogis(do.call("structure",c(list(.Data=unname(grpP[,1] + 2 * grpP[,2])), attributes(rp$GROUP_logitprobability))))
    dimnames(fit_P) <- dimnames(CIlow_P) <- dimnames(CIhigh_P) <- list(f$muNames[[3]],seq_len(nrow(data)))

    list(fit = fit_P,
         CI_low = CIlow_P,
         CI_high = CIhigh_P)
                    
}



##' Calculate group means for a data frame of covariates
##'
##' @param f 
##' @param data 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @export
getGroupMeans <- function(f, data, ...){
    UseMethod("getGroupMeans")
}


##' @rdname getGroupMeans
##' @method getGroupMeans mlld
##' @export
getGroupMeans.mlld <- function(f, data, keep.cov = FALSE, ...){
    require(lme4)
    if(missing(data))
        data <- f$data
    tmb_data <- f$tmb_data
    tmb_par <- f$pl
    tmb_map <- f$tmb_map
    tmb_random <- f$tmb_random
    tmb_profile <- f$tmb_profile

    tmb_data$X_pred <- Matrix::sparse.model.matrix(f$terms,
                                                   data = data,
                                                   transpose = TRUE,
                                                   row.names = FALSE,
                                                   xlev = f$xlevels)
    if(!inherits(tmb_data$X_pred,"TsparseMatrix"))
        tmb_data$X_pred <- as(tmb_data$X_pred,"TsparseMatrix")    

    tmb_data$XCom_pred <- Matrix::sparse.model.matrix(f$termsCommon,
                                                      data = data,
                                                      transpose = TRUE,
                                                      row.names = FALSE,
                                                      xlev = f$xlevelsCommon)
    if(!inherits(tmb_data$XCom_pred,"TsparseMatrix"))
        tmb_data$XCom_pred <- as(tmb_data$XCom_pred,"TsparseMatrix")    

    obj <- TMB::MakeADFun(data = tmb_data,
                     parameters = tmb_par,
                     map = tmb_map,
                     random = tmb_random,
                     profile = tmb_profile,
                     silent = TRUE,
                     DLL = "otoclass")
    rp <- obj$report()
    sdr <- TMB::sdreport(obj, obj$par, f$opt$he)
    ssdr <- TMB::summary.sdreport(sdr)    
    ##
    grpMean <- ssdr[rownames(ssdr) == "GROUP_mean_group",]
    fit_mean <- do.call("structure",c(list(.Data=unname(grpMean[,1])), attributes(rp$GROUP_mean_group)))
    dimnames(fit_mean) <- list(f$muNames[[2]],seq_len(nrow(data)),f$muNames[[3]])
    sd_mean <- do.call("structure",c(list(.Data=unname(grpMean[,2])), attributes(rp$GROUP_mean_group)))
    dimnames(sd_mean) <- list(f$muNames[[2]],seq_len(nrow(data)),f$muNames[[3]])
     ##
    grpCommon <- ssdr[rownames(ssdr) == "GROUP_mean_common",]
    fit_com <- do.call("structure",c(list(.Data=unname(grpCommon[,1])), attributes(rp$GROUP_mean_common)))
    dimnames(fit_com) <- list(f$muNames[[2]],seq_len(nrow(data)))
    sd_com <- do.call("structure",c(list(.Data=unname(grpCommon[,2])), attributes(rp$GROUP_mean_common)))
    dimnames(sd_com) <- list(f$muNames[[2]],seq_len(nrow(data)))
    ##
    grpTotal <- ssdr[rownames(ssdr) == "GROUP_mean_total",]
    fit_total <- do.call("structure",c(list(.Data=unname(grpTotal[,1])), attributes(rp$GROUP_mean_total)))
    dimnames(fit_total) <- list(f$muNames[[2]],seq_len(nrow(data)),f$muNames[[3]])
    sd_total <- do.call("structure",c(list(.Data=unname(grpTotal[,2])), attributes(rp$GROUP_mean_total)))
    dimnames(sd_total) <- list(f$muNames[[2]],seq_len(nrow(data)),f$muNames[[3]])
    ##
    r <- list(fit = list(group = fit_mean,
                         common = fit_com,
                         total = fit_total),
              se.fit = list(group = sd_mean,
                            common = sd_com,
                            total = sd_total))
    if(keep.cov){
        ii <- do.call("structure",c(list(.Data=unname(which(names(sdr$value) == "GROUP_mean_total"))), attributes(rp$GROUP_mean_total)))
        cov_list <- lapply(seq_len(dim(ii)[3]), function(gg){
            r <- lapply(seq_len(dim(ii)[2]), function(kk){
                sdr$cov[ii[,kk,gg],ii[,kk,gg]]
            })
            names(r) <- seq_len(nrow(data))
            r
        })
        names(cov_list) <- f$muNames[[3]]
        r$cov <- cov_list
    ##     cov_total <- sdr$cov[]
    ## rownames() <- colnames() <- 

    }
    return(r)
}


getGroupDataFrame <- function(x, ...){
    UseMethod("getGroupDataFrame")
}
getGroupDataFrame.mlld <- function(x,...){    
    as.data.frame(lapply(seq_along(x$groupLevels), function(i) factor(x$groupLevels[[i]][x$group+1],x$groupLevels[[i]])))
}

logspace_add <- function(logx,logy){
    pmax(logx, logy) + log1p(exp(-abs(logx - logy)))
}
logspace_sum <- function(logxs) Reduce(logspace_add,logxs)

rePrior <- function(x, prior){
    if(missing(prior) || is.null(prior))
        return(x)
    v <- x + log(prior)
    v - apply(v,2,logspace_sum)[col(v)]
}

##' @export
predict.mlld <- function(x, y = NULL, data = NULL, group = NULL, proportionGroup = NULL, prior = NULL, se.fit = FALSE, ...){
    pname <- ifelse(is.null(group), "posterior_logprobability_shape", "posterior_logprobability_all")
    if(!is.null(prior)){
        pname <- "posterior_logprobability_np"
    }
    if(is.null(y)){ ## Use fitted data
        if(!is.null(data))
            message("data is ignored when y is not given.")
        rp <- x$rp
        rn<- rownames(x$y)
        posterior <- t(exp(rePrior(rp[[pname]],prior)))
        fit <- aperm(rp$posterior_mean,c(2,1,3))
 
        if(se.fit){
            if(is.null(x$sdr))
                stop("The object must be fitted with doSdreport = TRUE to use se.fit")
            ssdr <- TMB::summary.sdreport(x$sdr)
        }

    }else{
        pname <- sprintf("pred_%s",pname)
        if(!is.matrix(y))
            if(is.vector(y)){
                y <- matrix(y,ncol = 1)
            }else{
                stop("y must be a matrix of observations")
            }        
        ## Check size of y:
        
        dat <- x$tmb_data
        if(is.null(data))
            data <- data.frame(ID = 1:nrow(y))
        dat$Y_pred <- t(y)
        if(length(setdiff(unique(c(x$all_vars,x$all_varsCommon,x$all_varsLogScale,x$all_varsTheta)), colnames(data))) > 0)
            warning("The new data is missing",)
##### Prepare model matrix #####
        dat$X_pred <- Matrix::sparse.model.matrix(x$terms,
                                                  data = model.frame(x$terms, data, na.action = na.pass),
                                                  transpose = TRUE,
                                                  row.names = FALSE,
                                                  xlev = x$xlevels)
        dat$X_pred <- as(dat$X_pred,"TsparseMatrix")

##### Prepare scale model matrix #####
        dat$XCom_pred <- Matrix::sparse.model.matrix(x$termsCommon,
                                                     data = model.frame(x$termsCommon,data, na.action = na.pass),
                                                     transpose = TRUE,
                                                     row.names = FALSE,
                                                     xlev = x$xlevelsCommon)
        dat$XCom_pred <- as(dat$XCom_pred,"TsparseMatrix")

        dat$XLogScale_pred <- Matrix::sparse.model.matrix(x$termsLogScale,
                                                     data = model.frame(x$termsLogScale,data, na.action = na.pass),
                                                     transpose = TRUE,
                                                     row.names = FALSE,
                                                     xlev = x$xlevelsLogScale)
        dat$XLogScale_pred <- as(dat$XLogScale_pred,"TsparseMatrix")

 ##### Prepare proportion model matrix #####
       dat$XTheta_pred <- Matrix::sparse.model.matrix(x$termsTheta,
                                                       data = model.frame(x$termsTheta,data, na.action = na.pass),
                                                       transpose = TRUE,
                                                       row.names = FALSE,
                                                       xlev = x$xlevelsTheta)
        dat$XTheta_pred <- as(dat$XTheta_pred,"TsparseMatrix")

        ## HANDLE RANDOM EFFECTS!
        
        ## if(is.null(proportionGroup)){# & is.null(prior)){
        ##     proportionGroup_pred <- rep(levels(x$proportionGroup)[1],ncol(dat$Y_pred))
        ##     ## }else if(is.null(proportionGroup) & !is.null(prior)){
        ##     ##     ## Change this to use a specified prior!!!
        ##     ##     proportionGroup_pred <- rep(levels(x$proportionGroup)[1],ncol(dat$Y_pred))
        ## }else{
        ##     proportionGroup_pred <- proportionGroup
        ## }
        ## proportionGroup_pred <- factor(proportionGroup_pred, levels = levels(x$proportionGroup))
        ## if(any(is.na(proportionGroup_pred)))
        ##     stop("proportionGroup must match the levels from the fitted object")

        ## pgp <- as.integer(proportionGroup_pred)-1
        ## attr(pgp,"levels") <- levels(x$proportionGroup)
        ## dat$proportionGroup_pred <- pgp

        if(is.null(group)){
            dat$G_pred <- matrix(NA_real_, ncol(dat$Y_pred), ncol(getGroupDataFrame(x)))
        }else{
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
            dat$G_pred <- group
            if(all(is.na(group)) || !isTRUE(all.equal(groupLevels, x$groupLevels)))
                warning("group does not match the fitted group.")
        }
        obj <- TMB::MakeADFun(dat,x$pl,x$tmb_map,
                              silent = x$silent, random = x$tmb_random,
                              DLL = "otoclass")
        rp <- obj$report(obj$env$last.par)
        rn <- names(y)

        posterior <- t(exp(rePrior(rp[[pname]],prior)))
        fit <- aperm(rp$pred_posterior_mean,c(2,1,3))

        if(se.fit){
            if(!is.null(x$sdr)){
                he <- solve(x$sdr$cov.fixed)
            }else{
                he <- NULL
            }
            sdr <- TMB::sdreport(obj,x$opt$par,he, getReportCovariance = FALSE)
            ssdr <- TMB::summary.sdreport(sdr)
        }
    }

    colnames(posterior) <- x$muNames[[3]]
    rownames(posterior) <- rn
    clas <- factor(colnames(posterior)[apply(posterior,1,which.max)], levels = x$muNames[[3]])
    dimnames(fit) <- list(rn,x$muNames[[2]],x$muNames[[3]])
    res <- list(class=clas,
                posterior = posterior,
                fit = fit)
    if(se.fit){
        if(is.null(x$sdr))
            stop("The object must be fitted with doSdreport = TRUE to use se.fit")
        se.fit <- aperm(array(ssdr[rownames(ssdr) == "posterior_mean",2],dim = c(dim(fit)[c(2,1,3)])),c(2,1,3))
        res$se.fit <- se.fit
    }
    
    return(res)
}


## predict.mlld <- function(x, type = c("class","mean"), newdata = NULL, dataTest = NULL, se.fit = FALSE, ...){
##     type <- match.arg(type)
##     if(type == "class"){
##         if(se.fit)
##             warning("se.fit is ignored for class predicitons")
##         if(is.null(newdata)){
##             return(list(class = x$test_class,
##                         posterior = x$test_posterior))
##         }else{
##             if(!is.matrix(newdata) & is.numeric(newdata) & ncol(newdata) == x$X_cols)
##                 stop(sprintf("newdata must be a numeric matrix with %d columns.",x$X_cols))
##             tt <- terms(x)
##             if(!isTRUE(all.equal(formula(terms(x)), ~1)) & is.null(dataTest))
##                 stop(sprintf("dataTest is required with columns: %s.", paste(x$all_vars,sep=", ")))
##             dat <- x$tmb_data
##             if(is.null(dataTest))
##                 dataTest <- data.frame(ID = 1:nrow(newdata))
##             dat$X_pred <- t(newdata[,x$features,drop=FALSE])
##             dat$covar_pred <- Matrix::sparse.model.matrix(tt,
##                                                           data = dataTest,
##                                                           transpose = TRUE,
##                                                           row.names = FALSE)
##             if(inherits(dat$covar_pred,"dgCMatrix"))
##                 dat$covar_pred <- as(dat$covar_pred,"dgTMatrix")
##             dat$estimateUnbiasedTestProportions <- 0
##             obj <- TMB::MakeADFun(dat,x$pl,x$tmb_map,
##                                   silent = x$silent, random = x$tmb_random,
##                                   DLL = "otoclass")
##             rp <- obj$report(obj$env$last.par)

##             posterior <- exp(t(rp$logpred))
##             colnames(posterior) <- x$muNames[[3]]
##             rownames(posterior) <- rownames(newdata)
##             clas <- factor(x$muNames[[3]][apply(posterior,1,which.max)], levels = x$muNames[[3]])
##             return(list(class = clas,
##                         posterior = posterior))
##         }
##     }else if(type == "mean"){
##         if(is.null(dataTest))
##             dataTest <- x$data
##         tt <- terms(x)
##         Xpred <- model.matrix(tt,
##                               data = dataTest,
##                               transpose = FALSE,
##                               row.names = FALSE)
##         mu <- x$rp$muUse
##         ttc <- x$termsCommon
##         Xpredc <- model.matrix(ttc,
##                               data = dataTest,
##                               transpose = FALSE,
##                               row.names = FALSE)
##         muc <- x$pl$commonMu
##         estMean <- sapply(1:dim(mu)[3], function(i){
##             v <- Xpred %*% mu[,,i] + Xpredc %*% muc
##             colnames(v) <- x$muNames[[2]]
##             rownames(v) <- rownames(dataTest)
##             return(v)
##         }, simplify = FALSE)
##         names(estMean) <- x$muNames[[3]]
##         if(!se.fit)
##             return(estMean)
##         if(is.null(x$sdr))
##             stop("The model must be fitted using doSdreport = TRUE")
##         indx <- which(names(x$sdr$value) %in% c("muUse","commonMu"))
##         names(indx) <- names(x$sdr$value)[indx]
##         isMu <- names(indx) == "mu"
##         covUse <- x$sdr$cov[indx,indx]
##         transMat <- cbind(Xpred,Xpredc)               
##         estSd <- sapply(1:dim(mu)[3], function(i){
##             v <- matrix(NA,nrow(Xpred),dim(mu)[2])
##             for(j in 1:dim(mu)[2]){
##                 indxIJ <- array(FALSE, dim(mu))
##                 indxIJ[,j,i] <- TRUE
##                 indxIJ <- as.vector(indxIJ)
##                 indxIJc <- matrix(FALSE, nrow(muc), ncol(muc))
##                 indxIJc[,j] <- TRUE
##                 v[,j] <- sqrt(diag(transMat %*% covUse[c(indxIJ,indxIJc),c(indxIJ,indxIJc)] %*% t(transMat)))
##             }
##             colnames(v) <- x$muNames[[2]]
##             rownames(v) <- rownames(dataTest)
##             return(v)
##         }, simplify = FALSE)
##         names(estSd) <- x$muNames[[3]]
##         return(list(fit = estMean, se.fit = estSd))
##     }
## }
