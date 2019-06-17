##' @importFrom stats predict vcov pnorm pt qnorm
##' @export
residuals.mlld <- function(object,...){
    pr <- stats::predict(object)
    Sigma <- stats::vcov(object)
    CC <- lapply(Sigma, function(ss) solve(t(chol(ss))))
    theta <- object$rp$theta
    Mvec <- object$rp$Mvec
    F <- sapply(1:nrow(object$y),function(i){
        th <- theta[,object$proportionGroup[i]]
        Muse <- Mvec[[object$confusionGroup[i]]][ifelse(is.na(object$group[i]),1,object$group[i]),]
        NC <- sum(Muse * th)
        Reduce("+",sapply(1:ncol(pr$posterior),function(j){
            y2 <- CC[[j]] %*% (object$y[i,,drop=TRUE] - pr$fit[i,,j])
            pv <- pkwmixE0(y2, object$rp$tmixp[,j], object$rp$df[,j],
                         object$rp$aParam[,j], object$rp$bParam[,j])
            ## pv <- (stats::pnorm(y2,0,1) * (1-object$rp$tmixp[,j]) + object$rp$tmixp[,j]*stats::pt(y2,object$rp$df[,j]))
            pv * th[j] * Muse[j] / NC            
        }, simplify = FALSE))
    })
    if(ncol(object$y) == 1)
        F <- matrix(F,nrow = 1)
    rownames(F) <- object$muNames[[2]]
    colnames(F) <- rownames(object$y)
    stats::qnorm(t(F))
}
