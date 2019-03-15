

##' @export
predict.mlld <- function(x, type = c("class","mean"), newdata = NULL, dataTest = NULL, se.fit = FALSE, ...){
    type <- match.arg(type)
    if(type == "class"){
        if(se.fit)
            warning("se.fit is ignored for class predicitons")
        if(is.null(newdata)){
            return(list(class = x$test_class,
                        posterior = x$test_posterior))
        }else{
            if(!is.matrix(newdata) & is.numeric(newdata) & ncol(newdata) == x$X_cols)
                stop(sprintf("newdata must be a numeric matrix with %d columns.",x$X_cols))
            tt <- terms(x)
            if(!isTRUE(all.equal(formula(terms(x)), ~1)) & is.null(dataTest))
                stop(sprintf("dataTest is required with columns: %s.", paste(x$all_vars,sep=", ")))
            dat <- x$tmb_data
            if(is.null(dataTest))
                dataTest <- data.frame(ID = 1:nrow(newdata))
            dat$X_pred <- t(newdata[,x$features,drop=FALSE])
            dat$covar_pred <- Matrix::sparse.model.matrix(tt,
                                                          data = dataTest,
                                                          transpose = TRUE,
                                                          row.names = FALSE)
            if(inherits(dat$covar_pred,"dgCMatrix"))
                dat$covar_pred <- as(dat$covar_pred,"dgTMatrix")
            dat$estimateUnbiasedTestProportions <- 0
            obj <- TMB::MakeADFun(dat,x$pl,x$tmb_map,
                                  silent = x$silent, random = x$tmb_random,
                                  DLL = "otoclass")
            rp <- obj$report(obj$env$last.par)

            posterior <- exp(t(rp$logpred))
            colnames(posterior) <- x$muNames[[3]]
            rownames(posterior) <- rownames(newdata)
            clas <- factor(x$muNames[[3]][apply(posterior,1,which.max)], levels = x$muNames[[3]])
            return(list(class = clas,
                        posterior = posterior))
        }
    }else if(type == "mean"){
        if(is.null(dataTest))
            dataTest <- x$data
        tt <- terms(x)
        Xpred <- model.matrix(tt,
                              data = dataTest,
                              transpose = FALSE,
                              row.names = FALSE)
        mu <- x$rp$muUse
        ttc <- x$termsCommon
        Xpredc <- model.matrix(ttc,
                              data = dataTest,
                              transpose = FALSE,
                              row.names = FALSE)
        muc <- x$pl$commonMu
        estMean <- sapply(1:dim(mu)[3], function(i){
            v <- Xpred %*% mu[,,i] + Xpredc %*% muc
            colnames(v) <- x$muNames[[2]]
            rownames(v) <- rownames(dataTest)
            return(v)
        }, simplify = FALSE)
        names(estMean) <- x$muNames[[3]]
        if(!se.fit)
            return(estMean)
        if(is.null(x$sdr))
            stop("The model must be fitted using doSdreport = TRUE")
        indx <- which(names(x$sdr$value) %in% c("muUse","commonMu"))
        names(indx) <- names(x$sdr$value)[indx]
        isMu <- names(indx) == "mu"
        covUse <- x$sdr$cov[indx,indx]
        transMat <- cbind(Xpred,Xpredc)               
        estSd <- sapply(1:dim(mu)[3], function(i){
            v <- matrix(NA,nrow(Xpred),dim(mu)[2])
            for(j in 1:dim(mu)[2]){
                indxIJ <- array(FALSE, dim(mu))
                indxIJ[,j,i] <- TRUE
                indxIJ <- as.vector(indxIJ)
                indxIJc <- matrix(FALSE, nrow(muc), ncol(muc))
                indxIJc[,j] <- TRUE
                v[,j] <- sqrt(diag(transMat %*% covUse[c(indxIJ,indxIJc),c(indxIJ,indxIJc)] %*% t(transMat)))
            }
            colnames(v) <- x$muNames[[2]]
            rownames(v) <- rownames(dataTest)
            return(v)
        }, simplify = FALSE)
        names(estSd) <- x$muNames[[3]]
        return(list(fit = estMean, se.fit = estSd))
    }
}
