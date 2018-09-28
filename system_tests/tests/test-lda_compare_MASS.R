##################################################################
######### Compare LDA calculations to the R package MASS #########
##################################################################


## Load packages
library(MASS)

set.seed(1)

## Test EFD calculation
invisible(
    replicate(20,{
        nclass <- sample(2:4,1)
        mu <- rnorm(nclass,0,3)
        sigma <- exp(rnorm(nclass,-1,0.1))
        S <- sample(1:length(mu),10000, replace = TRUE)
        X <- matrix(rnorm(length(S),mu[S],sigma[S]), ncol=1)
        ld1 <- lda(X,S)
        ld2 <- mlld(X,S, silent = TRUE)
        is_equal(as.vector(ld1$means),
                 ld2$opt$par[names(ld2$opt$par)=="mu"],
                 tolerance = 1e-4)
        p1 <- predict(ld1)
        is_equal(p1$posterior,
                 t(exp(ld2$rp$logpred)),
                 tolerance = 1e-3)
        is_equal(as.numeric(p1$class),
                 apply(ld2$rp$logpred,2,which.max))
    })
)

