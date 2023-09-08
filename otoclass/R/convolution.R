
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

medianpoolOne <- Vectorize(function(i,j,pic,rd,cd){
    mrow <- mean(1:rd)
    mcol <- mean(1:cd)
    im <- 1:rd - mrow
    jm <- 1:cd - mcol
    iu <- sapply(im,function(ii)min(max(1,i-ii),nrow(pic)))
    ju <- sapply(jm,function(ii)min(max(1,j-ii),ncol(pic)))
    median(pic[iu,ju])    
},c("i","j"))

medianpool <- function(pic,rd,cd){
    res <- matrix(NA,floor(nrow(pic)/rd),floor(ncol(pic)/cd))
    res <- outer(seq(1,nrow(pic),rd),
                 seq(1,ncol(pic),cd),
                 medianpoolOne,
                 pic = pic,
                 rd=rd,
                 cd=cd)
    return(res)
}

medianfilter <- function(pic, n){
    res <- outer(seq(1,nrow(pic),1),
                 seq(1,ncol(pic),1),
                 medianpoolOne,
                 pic = pic,
                 rd=n,
                 cd=n)
    return(res)
}

to01 <- function(x){
    x1 <- x - min(x)
    x1/max(x1)
}


edge_detect <- function(pic){
    k2 <- matrix(c(-1,-1,-1,-1,8,-1,-1,-1,-1),3,3)
    r <- convol(pic,k2)
    return(round(to01(r) * 255))
}

edge_detect2 <- function(pic){
    k2 <- matrix(c(1,0,-1,0,0,0,-1,0,1),3,3)
    r <- convol(pic,k2)
    return(round(to01(r) * 255))
}

edge_detect3 <- function(pic){
    k2 <- matrix(c(0,1,0,1,-4,1,0,1,0),3,3)
    r <- convol(pic,k2)
    return(round(to01(r) * 255))
}


sobel_x <- function(pic){    
    k2 <- matrix(c(-1,-2,-1,0,0,0,1,2,1),3,3)
    r <- convol(pic,k2)
    return(r)
}

sobel_y <- function(pic){    
    k2 <- t(matrix(c(-1,-2,-1,0,0,0,1,2,1),3,3))
    r <- convol(pic,k2)
    return(r)
}

sobel <- function(pic){
    r1 <- sobel_x(pic)
    r2 <- sobel_y(pic)
    r <- sqrt(r1^2 + r2^2)
    return(r)
}


sobel5_x <- function(pic){    
    k2 <- matrix(c(-5,-8,-10,-8,-5,-4,-10,-20,-10,-4,0,0,0,0,0,4,10,20,10,4,5,8,10,8,5),5,5)
    r <- convol(pic,k2)
    return(r)
}

sobel5_y <- function(pic){    
    k2 <- t(matrix(c(-5,-8,-10,-8,-5,-4,-10,-20,-10,-4,0,0,0,0,0,4,10,20,10,4,5,8,10,8,5),5,5))
    r <- convol(pic,k2)
    return(r)
}

sobel5 <- function(pic){
    r1 <- sobel5_x(pic)
    r2 <- sobel5_y(pic)
    r <- sqrt(r1^2 + r2^2)
    return(r)
}



prewitt_x <- function(pic){    
    k2 <- matrix(c(-1,-1,-1,0,0,0,1,1,1),3,3)
    r <- convol(pic,k2)
    return(r)
}

prewitt_y <- function(pic){    
    k2 <- t(matrix(c(-1,-1,-1,0,0,0,1,1,1),3,3))
    r <- convol(pic,k2)
    return(r)
}

prewitt <- function(pic){
    r1 <- prewitt_x(pic)
    r2 <- prewitt_y(pic)
    r <- sqrt(r1^2 + r2^2)
    return(r)
}


scharr_x <- function(pic){    
    k2 <- matrix(c(47,162,47,0,0,0,-47,-162,-47),3,3)
    r <- convol(pic,k2)
    return(r)
}

scharr_y <- function(pic){    
    k2 <- matrix(c(47,162,47,0,0,0,-47,-162,-47),3,3)
    r <- convol(pic,k2)
    return(r)
}

scharr <- function(pic){
    r1 <- scharr_x(pic)
    r2 <- scharr_y(pic)
    r <- sqrt(r1^2 + r2^2)
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

laplacian <- function(pic){
    k2 <- matrix(c(0,1,0,1,4,1,0,1,0),3,3)
    r <- convol(pic,k2)
    return(r)
}

sharpen <- function(pic){
    k2 <- matrix(c(0,-1,0,-1,5,-1,0,-1,0),3,3)
    r <- convol(pic,k2)
    return(r)
}

highpass <- function(pic){
    k2 <- matrix(c(0,-1/4,0,-1/4,2,-1/4,0,-1/4,0),3,3)
    r <- convol(pic,k2)
    return(r)
}

meanfilter <- function(pic, n = 3){
    k2 <- matrix(1/n^2,n,n)
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

between <- function(x,a,b,leq=TRUE,req=TRUE){
    a < x && x < b || (leq && a==x) || (req && x==b)
}

canny_edge <- function(pic){
    ## Blur
    p <- gaussian_blur(extendPic(pic,21), 21, 15)
    ## Sobel
    Ix <- sobel5_x(extendPic(p,5))
    Iy <- sobel5_y(extendPic(p,5))
    G <- sqrt(Ix^2 + Iy^2)
    Th <- atan2(Iy,Ix)
    ## Non-maximum suppression
    angle <- Th * 180 / pi
    angle[angle < 0] <- angle[angle < 0] + 180
    getOne <- Vectorize(function(i,j){
        q <- 255
        r <- 255
        ##angle 0
        if( between(angle[i,j],0,22.5,TRUE,FALSE) || between(angle[i,j],157.5,180,TRUE,TRUE) ){
            q <- G[i, pmin(j+1,ncol(G))]
            r <- G[i, pmax(j-1,1)]
            ##angle 45
        }else if(between(angle[i,j],22.5,67.5,TRUE,FALSE)){
            q <- G[pmin(i+1,nrow(G)), pmax(j-1,1)]
            r <- G[pmax(i-1,1), pmin(j+1,ncol(G))]
            ##angle 90
        }else if(between(angle[i,j],67.5,112.5,TRUE,FALSE)){
            q <- G[pmin(i+1,nrow(G)), j]
            r <- G[pmax(i-1,1), j]
            ##angle 135
        }else if(between(angle[i,j],112.5,157.5,TRUE,FALSE)){
            q <- G[pmax(i-1,1), max(j-1,1)]
            r <- G[pmin(i+1,nrow(G)), pmin(j+1,ncol(G))]
        }else{
            q <- r <- 255
        }
        if( (G[i,j] >= q) && (G[i,j] >= r)){
            return(G[i,j])
        }
        return(0)
    })
    p2 <- outer(seq_len(nrow(pic)),seq_len(ncol(pic)), getOne)
    ## Double treshold
    lowThresholdRatio <- 0.05
    highThresholdRatio <- 0.09
    highThreshold <- 0.3 # max(p2) * highThresholdRatio;
    lowThreshold <- 0.1 # highThreshold * lowThresholdRatio;
    p3 <- 0 * p2
    p3[p2 >= lowThreshold] <- 25
    p3[p2 >= highThreshold] <- 255
    ## Hysteresis
    doOneHy <- Vectorize(function(i,j){
        if(i == 0 || j == 0 || i == nrow(p3) || j == ncol(p3))
            return(0)
        if(p3[i,j] != 25)
            return(p3[i,j])
        if(any(p3[(i-1):(i+1), (j-1):(j+1)]==255))
            return(255)
        return(0)
    })
    p4 <- outer(seq_len(nrow(pic)),seq_len(ncol(pic)), doOneHy)
    p4
}
