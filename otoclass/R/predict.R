

##' @export
predict.mlld <- function(x, type = c("class","mean"), newdata = NULL, dataTest = NULL, ...){
    type <- match.arg(type)
    if(type == "class"){
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
        mu <- x$pl$mu
        estMean <- sapply(1:dim(mu)[3], function(i){
            v <- Xpred %*% mu[,,i]
            colnames(v) <- x$muNames[[2]]
            rownames(v) <- rownames(dataTest)
            return(v)
        }, simplify = FALSE)
        names(estMean) <- x$muNames[[3]]
        return(estMean)
    }
}
