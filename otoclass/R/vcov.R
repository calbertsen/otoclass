##' @export
vcov.mlld <- function(object,...){
    ## sigmas <- exp(x$pl$logSigma)
    ## corpars <- x$pl$corpar
    ## Sigma <- sapply(1:ncol(sigmas), function(i) {
    ##     A <- diag(sigmas[, i], nrow(sigmas), nrow(sigmas))
    ##     M <- diag(1, nrow(A), nrow(A))
    ##     M[lower.tri(M)] <- corpars[, i]
    ##     A %*% cov2cor(M %*% t(M)) %*% t(A)
    ## }, simplify = FALSE)
    Sigma <- lapply(object$rp$SigmaList,function(xx){
        yy <- xx
        rownames(yy) <- colnames(yy) <- object$muNames[[2]]
        yy
        })
    names(Sigma) <- object$muNames[[3]]
    Sigma
}

