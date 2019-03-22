
##' @export
vcov.mlld <- function(x,...){
    ## sigmas <- exp(x$pl$logSigma)
    ## corpars <- x$pl$corpar
    ## Sigma <- sapply(1:ncol(sigmas), function(i) {
    ##     A <- diag(sigmas[, i], nrow(sigmas), nrow(sigmas))
    ##     M <- diag(1, nrow(A), nrow(A))
    ##     M[lower.tri(M)] <- corpars[, i]
    ##     A %*% cov2cor(M %*% t(M)) %*% t(A)
    ## }, simplify = FALSE)
    Sigma <- lapply(x$rp$SigmaList,function(xx){
        yy <- xx
        rownames(yy) <- colnames(yy) <- x$muNames[[2]]
        yy
        })
    names(Sigma) <- x$muNames[[3]]
    Sigma
}


##' @export
residuals.mlld <- function(x,...){
    pr <- predict(x)
    Sigma <- vcov(x)
    CC <- lapply(Sigma, function(ss) solve(t(chol(ss))))
    theta <- x$rp$theta
    Mvec <- x$rp$Mvec
    F <- sapply(1:nrow(x$y),function(i){
        th <- theta[,x$proportionGroup[i]]
        Muse <- Mvec[[x$confusionGroup[i]]][ifelse(is.na(x$group[i]),1,x$group[i]),]
        NC <- sum(Muse * th)
        Reduce("+",sapply(1:ncol(pr$posterior),function(j){
            y2 <- CC[[j]] %*% (x$y[i,,drop=TRUE] - pr$fit[i,,j])
            (pnorm(y2,0,1) * (1-x$rp$tmixp[,j]) + x$rp$tmixp[,j]*pt(y2,x$rp$df[,j])) * th[j] * Muse[j] / NC            
        }, simplify = FALSE))
    })
    if(ncol(x$y) == 1)
        F <- matrix(F,nrow = 1)
    rownames(F) <- x$muNames[[2]]
    colnames(F) <- rownames(x$y)
    qnorm(t(F))
}
