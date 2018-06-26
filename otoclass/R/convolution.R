
## convOne <- Vectorize(function(i,j,pic,kernel){
##     mrow <- mean(1:nrow(kernel))
##     mcol <- mean(1:ncol(kernel))
##     im <- 1:nrow(kernel) - mrow
##     jm <- 1:ncol(kernel) - mcol
##     iu <- sapply(im,function(ii)min(max(1,i-ii),nrow(pic)))
##     ju <- sapply(jm,function(ii)min(max(1,j-ii),ncol(pic)))
##     sum(pic[iu,ju] * kernel)
##     ## sum(outer(1:nrow(kernel) - mrow,
##     ##       1:ncol(kernel) - mcol,
##     ##       Vectorize(function(im,jm){
##     ##           iu <- min(max(1,i-im),nrow(pic))
##     ##           ju <- min(max(1,j-jm),ncol(pic))
##     ##           pic[iu,ju] * kernel[im + mrow,jm + mcol]
##     ##       })))
##    },c("i","j"))

## conv <- function(pic,kernel)
##     outer(1:nrow(pic),
##           1:ncol(pic),
##           convOne,
##           pic = pic,
##           kernel = kernel)


convol <- function(pic, kernel){
    .Call("convol2d",x = pic, k = kernel, PACKAGE = "otoclass")
}

maxpoolOne <- Vectorize(function(i,j,pic,rd,cd){
    mrow <- mean(1:rd)
    mcol <- mean(1:cd)
    im <- 1:rd - mrow
    jm <- 1:cd - mcol
    iu <- sapply(im,function(ii)min(max(1,i-ii),nrow(pic)))
    ju <- sapply(jm,function(ii)min(max(1,j-ii),ncol(pic)))
    max(pic[iu,ju])    
},c("i","j"))

maxpool <- function(pic,rd,cd){
    res <- matrix(NA,floor(nrow(pic)/rd),floor(ncol(pic)/cd))
    res <- outer(seq(1,nrow(pic),rd),
                 seq(1,ncol(pic),cd),
                 maxpoolOne,
                 pic = pic,
                 rd=rd,
                 cd=cd)
    return(res)
}

to01 <- function(x){
    x1 <- x - min(x)
    x1/max(x1)
}


edge_detect <- function(pic){
    k2 <- matrix(c(-1,-1,-1,-1,8,-1,-1,-1,-1),3,3)
    r <- convol(pic,k2)
    return(r)
}

edge_detect2 <- function(pic){
    k2 <- matrix(c(1,0,-1,0,0,0,-1,0,1),3,3)
    r <- convol(pic,k2)
    return(r)
}

edge_detect3 <- function(pic){
    k2 <- matrix(c(0,1,0,1,-4,1,0,1,0),3,3)
    r <- convol(pic,k2)
    return(r)
}


sobel <- function(pic){    
    k2 <- matrix(c(-1,-2,-1,0,0,0,1,2,1),3,3)
    r <- convol(pic,k2)
    return(r)
}



##' @importFrom stats median dnorm
gaussian_blur <- function(pic,n = 3, sigma = 1){
    xx <- 1:n - stats::median(1:n)
    k2 <- outer(xx,xx, function(x,y) stats::dnorm(x,0,sigma) * stats::dnorm(y,0,sigma))
    r <- convol(pic,k2/sum(k2))
    return(r)
}

##' @importFrom stats median dnorm
unsharpmask <- function(pic, n = 3, sigma = 1){
    ## k2 <- -1/256 * matrix(c(1,4,6,4,1,
    ##                         4,16,24,16,4,
    ##                         6,24,-476,24,6,
    ##                         4,16,24,16,4,
    ##                         1,4,6,4,1),
    ##                       5,5, byrow = TRUE)
    xx <- 1:n - stats::median(1:n)
    k2 <- outer(xx,xx, function(x,y) stats::dnorm(x,0,sigma) * stats::dnorm(y,0,sigma))
    s <- sum(k2)
    k2[stats::median(1:n),stats::median(1:n)] <- k2[stats::median(1:n),stats::median(1:n)] - 2 * s
    r <- convol(pic,-k2/s)
    return(r) ## round(plogis(r-mean(r)) * 255)
}

highpass <- function(pic){
    k2 <- matrix(c(0,-1/4,0,-1/4,2,-1/4,0,-1/4,0),3,3)
    r <- convol(pic,k2)
    return(r)
}

meanfilter <- function(pic){
    k2 <- matrix(1/9,3,3)
    r <- convol(pic,k2)
    return(r)
}

lowpass <- function(pic){
    k2 <- matrix(c(0,1/8,0,1/8,1/2,1/8,0,1/8,0),3,3)
    r <- convol(pic,k2)
    return(r)
}


emboss <- function(pic){
    k2 <- matrix(c(-2,-1,0,-1,1,1,0,1,2),3,3)
    r <- convol(pic,k2)
    return(r)
}
