

##' @export
predict.mlld <- function(x, y = NULL, data = NULL,proportionGroup = NULL, prior = NULL, se.fit = FALSE, ...){
    if(is.null(y)){ ## Use fitted data
        if(!is.null(data))
            message("data is ignored when y is not given.")
        rp <- x$rp
        rn<- rownames(x$y)
        posterior <- t(exp(rp$posterior_logprobability))
        fit <- aperm(rp$posterior_mean,c(2,1,3))
 
        if(se.fit){
            if(is.null(x$sdr))
                stop("The object must be fitted with doSdreport = TRUE to use se.fit")
            ssdr <- TMB::summary.sdreport(x$sdr)
        }

    }else{
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
        dat$X_pred <- Matrix::sparse.model.matrix(x$terms,
                                                  data = data,
                                                  transpose = TRUE,
                                                  row.names = FALSE,
                                                  xlev = x$xlevels)
        if(inherits(dat$X_pred,"dgCMatrix"))
            dat$X_pred <- as(dat$X_pred,"dgTMatrix")
        dat$XCom_pred <- Matrix::sparse.model.matrix(x$termsCommon,
                                                     data = data,
                                                     transpose = TRUE,
                                                     row.names = FALSE,
                                                     xlev = x$xlevelsCommon)
        if(inherits(dat$XCom_pred,"dgCMatrix"))
            dat$XCom_pred <- as(dat$XCom_pred,"dgTMatrix")
        if(is.null(proportionGroup) & is.null(prior)){
            proportionGroup_pred <- rep(levels(x$proportionGroup)[1],ncol(dat$Y_pred))
        }else if(is.null(proportionGroup) & !is.null(prior)){
            ## Change this to use a specified prior!!!
            proportionGroup_pred <- rep(levels(x$proportionGroup)[1],ncol(dat$Y_pred))
        }else{
            proportionGroup_pred <- proportionGroup
        }
        proportionGroup_pred <- factor(proportionGroup_pred, levels = levels(x$proportionGroup))
        if(any(is.na(proportionGroup_pred)))
            stop("proportionGroup must match the levels from the fitted object")

        pgp <- as.integer(proportionGroup_pred)-1
        attr(pgp,"levels") <- levels(x$proportionGroup)
        dat$proportionGroup_pred <- pgp

        obj <- TMB::MakeADFun(dat,x$pl,x$tmb_map,
                              silent = x$silent, random = x$tmb_random,
                              DLL = "otoclass")
        rp <- obj$report(obj$env$last.par)
        rn <- names(y)

        posterior <- t(exp(rp$pred_posterior_logprobability))
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
    clas <- colnames(posterior)[apply(posterior,1,which.max)]
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
