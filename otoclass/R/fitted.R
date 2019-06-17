##' @export
fitted.mlld <- function(object,...){
   pr <- predict(object)
    Sigma <- vcov(object)
    CC <- lapply(Sigma, function(ss) solve(t(chol(ss))))
    theta <- object$rp$theta
    Mvec <- object$rp$Mvec
    E <- sapply(1:nrow(object$y), function(i) {
        th <- theta[, object$proportionGroup[i]]
        Muse <- Mvec[[object$confusionGroup[i]]][ifelse(is.na(object$group[i]), 
            1, object$group[i]), ]
        NC <- sum(Muse * th)
        Reduce("+", sapply(1:ncol(pr$posterior), function(j) {
            pr$fit[i,, j] * th[j] * Muse[j]/NC
        }, simplify = FALSE))
    })
   if (ncol(object$y) == 1) 
       E <- matrix(E, nrow = 1)
   rownames(E) <- object$muNames[[2]]
   colnames(E) <- rownames(object$y)
   t(E)
}
