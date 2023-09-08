##' @export
fitted.mlld <- function(x,...){
    ## pr <- predict(x)
    ##  Sigma <- vcov(x)
    ##  CC <- lapply(Sigma, function(ss) solve(t(chol(ss))))
    ##  theta <- x$rp$theta
    ##  Mvec <- x$rp$Mvec
    ##  E <- sapply(1:nrow(x$y), function(i) {
    ##      th <- theta[, x$proportionGroup[i]]
    ##      Muse <- Mvec[[x$confusionGroup[i]]][ifelse(is.na(x$group[i]), 
    ##          1, x$group[i]), ]
    ##      NC <- sum(Muse * th)
    ##      Reduce("+", sapply(1:ncol(pr$posterior), function(j) {
    ##          pr$fit[i,, j] * th[j] * Muse[j]/NC
    ##      }, simplify = FALSE))
    ##  })
    ## if (ncol(x$y) == 1) 
    ##     E <- matrix(E, nrow = 1)
    E0 <- x$rp$posterior_mean
    P_S <- plogis(x$rp$prior_logitprobability)
    Px <- aperm(array(P_S,dim=c(nrow(P_S),ncol(P_S),dim(E0)[1])),c(3,2,1))
    E <- apply(E * Px,1:2,sum)
    rownames(E) <- x$muNames[[2]]
    colnames(E) <- rownames(x$y)
    t(E)
}
