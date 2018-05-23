library(TMB)
library(otoclass)
library(Matrix)

nobs <- 100
npred <- 20
nxtra <- 100
X1 <- matrix(c(rnorm(nobs,-1,2),
       rnorm(nobs,1,1)),ncol=1)
G <- rep(c("A","B"), each = nobs)

X <- cbind(X1,matrix(rnorm(nrow(X1)*nxtra,0,5),ncol=nxtra))

X2 <- matrix(c(rnorm(npred,-1,2),
       rnorm(npred,1,1)),ncol=1)
Gp <- rep(c("A","B"), each = npred)
Xp <- cbind(X2,matrix(rnorm(nrow(X2)*nxtra,0,5),ncol=nxtra))

Xscale <- scale(X,rep(0,ncol(X)),rep(1,ncol(X)))
dat <- list(model = 0,
            X = t(Xscale),
            G = factor(G),
            Q = as(matrix(0.0,0,0),"sparseMatrix"),
            penalty = 1,
            prior = c(0.5,0.5),
            X_pred = t(scale(Xp,attr(Xscale,"scaled:center"),attr(Xscale,"scaled:scale"))),
            logLambda = log(0.4)
            )
n <- nrow(dat$X)
par <- list(mu = matrix(apply(dat$X,1,mean),n,nlevels(dat$G)),
            efd = array(0.0,dim=c(0,0,0)),
            logSigma = matrix(log(apply(dat$X,1,sd)),n,nlevels(dat$G)),
            corpar = matrix(0,(n*n-n)/2,nlevels(dat$G)), ## t(chol(cor(X)))[lower.tri(t(chol(cor(X))))]
            logDelta = 0
           )
map <- list(logSigma = factor(row(par$logSigma)),
            corpar = factor(rep(NA,length(par$corpar))), # factor(row(par$corpar)) # factor(rep(NA,length(par$corpar))) #
            logDelta = factor(NA)
            )

obj <- MakeADFun(dat,par,map, DLL = "otoclass")

obj$fn()
obj$env$data$logLambda <- log(1)
obj$fn()
obj$env$data$logLambda <- log(0.01)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 100000, eval.max = 100000))
opt

tab <- table(Gp,apply(obj$report(opt$par)$logpred,2,which.max))
tab
sum(diag(tab)) / sum(tab)

penv <- round(exp(obj$env$parList(x=opt$par)$logSigma) * (obj$env$parList(x=opt$par)$mu - rowMeans(obj$env$parList(x=opt$par)$mu)),4)

plot(abs(penv[,1]))


doCV <- function(l){
    cat(l,"\n")
    obj$env$data$logLambda <- log(l)
    capture.output(opt2 <- nlminb(opt$par, obj$fn, obj$gr, control = list(iter.max = 100000, eval.max = 100000)))
    tab <- table(Gp,apply(obj$report(opt2$par)$logpred,2,which.max))
    sum(diag(tab)) / sum(tab)
}

cvvals <- sapply(seq(0.01,0.2,0.05),doCV)


library(penalizedLDA)

re <- PenalizedLDA((X),as.numeric(factor(G)),(Xp),K=1,lambda=0.4)

tab2 <- table(Gp,re$ypred)
tab2
sum(diag(tab2)) / sum(tab2)


library(MASS)

lr <- lda(X,G, c(0.5,0.5))
tab3 <- table(Gp,predict(lr,newdata=Xp)$class)
sum(diag(tab3)) / sum(tab3)



lr <- lda(X[,1,drop=FALSE],G, c(0.5,0.5))
tab4 <- table(Gp,predict(lr,newdata=Xp[,1,drop=FALSE])$class)
sum(diag(tab4)) / sum(tab4)


### FS
library(otoclass)

indx <- fcbf(X,factor(G),0.1)
od <- otoclass::discrim(X[,indx,drop=FALSE],factor(G),Xp[,indx,drop=FALSE])

tab5 <- table(Gp,od$predicted)
tab5
sum(diag(tab5)) / sum(tab5)


### KNN

library(otoclass)

k2 <- knn(X[,1:5],factor(G),Xp[,1:5],10,"L2")

tab6 <- table(Gp,k2$predicted)
tab6
sum(diag(tab6)) / sum(tab6)


###########################################################################
###########################################################################
################################ GMRF TEST ################################
###########################################################################
###########################################################################

library(TMB)
library(otoclass)
library(Matrix)

im <- read_image("~/dtu/otolith_image/Aku13A/J0328_00001_0_5x1x2_0x1_1.jpg",
                 logisticTransform = TRUE)
plot(im,col="red")
grp1 <- efd(im[[1]],60, normalize = TRUE, returnAsList = TRUE)
grp1$A0 <- grp1$C0 <- NULL

grp2 <- grp1
grp2$D[1] <- htp2$D[1] + 0.01
grp2$C[20] <- htp2$C[20] + 0.01

nOto <- 100

coordGrp1 <- do.call("efd2coord",c(n=nOto,grp1))
coordGrp2 <- do.call("efd2coord",c(n=nOto,grp2))
plot(coordGrp1, type = "l", asp = 1)
lines(coordGrp2, col = "red")



### Simulate otoliths
phi <- 0.0001
## Setup Covariance structure
cphi <- 1/phi+phi-2
pnts <- c(1,1,2+cphi)
t <- seq(0,1,len = nOto)
for(i in 2:length(t)){
    pnts <- rbind(pnts,
                  c(i,i,2+cphi),
                  c(i,i-1,-1),
                  c(i-1,i,-1)
                  )
}
pnts <- rbind(pnts,
              c(length(t),1,-1),
              c(1,length(t),-1)
              )
Q <- sparseMatrix(i=pnts[,1],j=pnts[,2],x=pnts[,3])
Qsim <- Q
diag(Qsim) <- diag(Qsim) + 0.5
sigma <- as.matrix(solve(Qsim))
sdx <- 0.005/sqrt(sigma[1,1])
sdy <- 0.025/sqrt(sigma[1,1])
rho2 <- 0.0
sigma2 <- matrix(c(sdx^2,sdx*sdy*rho2,sdx*sdy*rho2,sdy^2),ncol=2)
sepext <- Matrix::kronecker(sigma,sigma2)


nSim <- 120

X1 <- replicate(nSim,matrix(argosTrack:::rmvnorm(1,
                                               as.vector(t(coordGrp1)),
                                               sepext),nrow=2), simplify = FALSE)
X2 <- replicate(nSim,matrix(argosTrack:::rmvnorm(1,
                                  as.vector(t(coordGrp2)),
                                  sepext),nrow=2), simplify = FALSE)
X <- simplify2array(c(X1,X2))
G <- rep(1:2,each=nSim)


plot(coordGrp1, type = "l", asp = 1)
lines(coordGrp2, col = "red")
invisible(sapply(1:dim(X)[3],function(i)lines(t(X[,,i]), lty = 1, col = c(rgb(0,0,0,0.9),rgb(1,0,0,0.9))[G[i]])))

## Fit model
dat <- list(model = 2,
            X = X,
            G = factor(G),
            Q = Q,
            penalty = 2,
            prior = c(0.5,0.5),
            X_pred = X,
            logLambda = log(0.4)
            )
n <- nrow(dat$X)
efds <- array(0,dim = c(4,60,nlevels(dat$G)))
efds[] <- runif(prod(dim(efds)),-0.01,0.01)
efds[1,1,] <- 1
efds[2:3,1,] <- 0
efds[4,1,] <- 1
efdMap <- array(1:prod(dim(efds)),dim = dim(efds))
efdMap[1:3,1,] <- NA
par <- list(mu = matrix(0,0,0),
            efd = efds,
            logSigma = matrix(0,n,nlevels(dat$G)),
            corpar = matrix(0,(n*n-n)/2,nlevels(dat$G)),
            logDelta = 0
           )
map <- list(logSigma = factor(row(par$logSigma)),
            corpar = factor(rep(NA,length(par$corpar))),
            efd = factor(efdMap)
            )

obj <- MakeADFun(dat,par,map, DLL = "otoclass")


obj$fn()
obj$env$data$logLambda <- log(400)
obj$fn()
obj$env$data$logLambda <- log(.4)
obj$fn()


opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 100000, eval.max = 100000))
opt

pl <- obj$env$parList(opt$par)

round(pl$efd,4)

rp <- obj$report(opt$par)


plot(coordGrp1, type = "l", asp = 1)
lines(coordGrp2, col = "red")
sapply(1:2,function(i) lines(t(rp$coords[[i]]),col=i,lty = i+1,lwd=4))

round(obj$env$parList(opt$par)$efd,3)

round(obj$env$parList(opt$par)$efd[,,1] - obj$env$parList(opt$par)$efd[,,2],3)

table(G,apply(rp$logpred,2,which.max))


### Compare LDA

ef <- t(apply(aperm(X,c(2,1,3)),3,efd, N=60, normalize = TRUE))
noUse <- grep("(A1|B1|C1|A0|C0)",colnames(ef))

ldr <- MASS::lda(ef[,-noUse],G,c(0.5,0.5))
table(G,predict(ldr,ef[,-noUse])$class)



library(penalizedLDA)

re <- PenalizedLDA((ef[,-noUse]),as.numeric(factor(G)),ef[,-noUse],K=1,lambda=.1)

tab2 <- table(G,re$ypred)
tab2
sum(diag(tab2)) / sum(tab2)
