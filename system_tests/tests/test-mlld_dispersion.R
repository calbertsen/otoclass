##################################################################
################# Test dispersion group argument #################
##################################################################

set.seed(234)
library(otoclass)

## No dispersion
N <- 4000
Grp <- rbinom(N,1,0.3)
X <- rnorm(N, 0.2 + 0.5 * Grp, 1)

fit <- mlld(matrix(X,ncol=1),
            factor(Grp))

fit

fit$rp$theta
fit$rp$thetaDisp


dg <- factor(rep(1:4,len = N))
dg[1:20] <- NA_real_
fit <- mlld(matrix(X,ncol=1),
            factor(Grp),
            dispersionGroup=dg)

fit

fit$rp$theta

fit$rp$thetaDisp



## No dispersion
rm("fit")
N <- 10000
DG <- sample(1:10,N,TRUE)
d <- rnorm(max(DG),0,0.1)
Grp <- rbinom(N,1,plogis(qlogis(0.3) + d[DG]))

apply(sapply(split(Grp,DG),table),2,function(x)x/sum(x))

X <- rnorm(N, 0.2 + 0.5 * Grp, 1)

fit <- mlld(matrix(X,ncol=1),
            factor(Grp),
            dispersionGroup=factor(DG),
            doSdreport=TRUE,
            getReportCovariance=FALSE)
fit

fit$rp$thetaDisp

fit$rp$theta

ssdr <- summary(fit$sdr)

ssdr[rownames(ssdr) %in% "theta",]


fit2 <- mlld(matrix(X,ncol=1),
            factor(Grp),
            doSdreport=TRUE,
            getReportCovariance=FALSE)
fit2


ssdr2 <- summary(fit2$sdr)

ssdr2[rownames(ssdr2) %in% "theta",]
