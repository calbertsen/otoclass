
##' @importFrom stats approx
setNumObs_image<- function(dat,n){
    t <- seq(0,1,len=n)
    tmpX <- stats::approx(x=seq(0,1,len=dim(dat)[1]),y=dat[,1],xout=t)
    tmpY <- stats::approx(x=seq(0,1,len=dim(dat)[1]),y=dat[,2],xout=t)
    datOut <- cbind(tmpX$y,tmpY$y)
    attrs <- attributes(dat)
    attrs[names(attributes(datOut))] <- attributes(datOut)
    attributes(datOut) <- attrs
    return(datOut)
}
center_image <- function(dat){
    tmpX <- dat[,1]
    tmpY <- dat[,2]
    datOut <- cbind(tmpX-mean(range(tmpX)),
                    tmpY-mean(range(tmpY)))
    attrs <- attributes(dat)
    attrs[names(attributes(datOut))] <- attributes(datOut)
    attributes(datOut) <- attrs
    return(datOut)
}
scale_image <- function(dat){
    tmpX <- dat[,1]
    tmpY <- dat[,2]
    datOut <- cbind(tmpX/diff(range(tmpX)),
                 tmpY/diff(range(tmpX)))
    attrs <- attributes(dat)
    attrs[names(attributes(datOut))] <- attributes(datOut)
    attributes(datOut) <- attrs
    return(datOut)
}

findMaxDist_image <- function(dat){
   xsq <- outer(dat[,1],dat[,1],function(x,y)(x-y)^2)
   ysq <- outer(dat[,2],dat[,2],function(x,y)(x-y)^2)
   val <- sqrt(xsq+ysq)
   return(which(val==max(val),arr.ind=TRUE)[1,])
}

rotate_image <- function(dat){
    oldCenter <- apply(dat,2,mean)
    dat <- center_image(dat)
    newVecIndx <- findMaxDist_image(dat)
    newVec <- dat[newVecIndx[1],]-dat[newVecIndx[2],]
    if(newVec[2]<0)
        newVec <- -newVec
    alpha <- acos(newVec[1]/(sqrt(sum(newVec^2))))
    rotMat <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)), ncol=2)
    val <- dat%*%rotMat
    datOut <- t(apply(val,1,function(x)x+oldCenter))
    attrs <- attributes(dat)
    attrs[names(attributes(datOut))] <- attributes(datOut)
    attributes(datOut) <- attrs
    return(datOut)
}

newStart_image <- function(dat){
    ldat <- dim(dat)[1]
    newBegin <- which.max(dat[,1])
    if(newBegin == 1){
        datNew <- dat
    }else{
        datNew <- rbind(dat[newBegin:ldat,],
                        dat[1:(newBegin-1),]
                        )
    }
    if(datNew[2,2] < datNew[1,2])
        datNew <- rbind(datNew[1,],
                        apply(datNew[-1,],2,rev)
                        )
    attrs <- attributes(dat)
    attrs[names(attributes(datNew))] <- attributes(datNew)
    attributes(datNew) <- attrs
    return(datNew)
}

flip_image <- function(dat,datCompare,forceFlip=FALSE){
    # Both dat and datCompare must be preprocessed
    dat1 <- dat
    dat2 <- dat
    dat2[,2] <- -dat2[,2]
    dat2 <-  newStart_image(dat2)
    attr(dat2,"Flipped") <- TRUE

    if(forceFlip)
        return(dat2)
    
    dat1tmp <- setNumObs_image(dat1,1000)
    dat2tmp <- setNumObs_image(dat2,1000)
    datComptmp <- setNumObs_image(datCompare,1000)
    mse1 <- mean((dat1tmp-datComptmp)^2)
    mse2 <- mean((dat2tmp-datComptmp)^2)
   
    if(mse1 <= mse2){
        return(dat1)
    }else{
        return(dat2)
    }
}

LineLength <- function(v1,v2){
    sqrt(sum((v1-v2)^2))
}
##' Get a pixel matrix from an image file
##'
##' @param file Path to image file
##' @param grey Should output be greyscale?
##' @return A matrix of pixel values (0-255)
##' @author Christoffer Moesgaard Albertsen
##' @importFrom jpeg readJPEG
##' @importFrom png readPNG
##' @importFrom tiff readTIFF
getPixelMatrix <- function(file, grey=TRUE){
    if(grepl("\\.jp(e)*g$",file,ignore.case=TRUE)){
        im <- jpeg::readJPEG(file)
    }else if(grepl("\\.png$",file,ignore.case=TRUE)){
        im <- png::readPNG(file)
    }else if(grepl("\\.tif(f)*$",file,ignore.case=TRUE)){
        im <- tiff::readTIFF(file)
    }else{
        stop(sprintf("Unsupported file format: %s. Only JPEG, PNG, and TIFF are supported.",tail(strsplit(file,".")[[1]],1)))
    }
    if(grey & length(dim(im)) == 3){
        if(dim(im)[3] >= 3)
            imOut <- (0.2989 * im[,,1] + 0.5870 * im[,,2] + 0.1140 * im[,,3]) * 255
    }else{
        imOut <- im * 255
    }
    return(imOut)
}

##' Read Otolith Images and Extract Contours
##'
##' @param file 
##' @param noiseFactor 
##' @param onlyOne 
##' @param minPixelDiff 
##' @param extreme 
##' @param pixelwise 
##' @param assignSinglesByPosition 
##' @return otolith image information
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats kmeans relevel
##' @export
read_image<- function(file,noiseFactor = NULL, onlyOne = FALSE, minPixelDiff = 0.05 * min(nc,nr), extreme = TRUE, pixelwise = FALSE, assignSinglesByPosition = TRUE){
    r<-getPixelMatrix(file)
    rv <- t(r[nrow(r):1,])
    maxrv <- max(rv)
    nc <- ncol(rv)
    nr <- nrow(rv)

    whiteBorder <- mean(c(rv[c(1:(nr*0.02),nr:(nr-nr*0.02)),],rv[,c(1:(nc*0.02),nc:(nc-nc*0.02))])) > 255/2
    if(whiteBorder){
        rv <- 255 - rv
    }

    if(is.null(noiseFactor)){
        hh<-hist(rv,breaks=(-1):maxrv + 0.5,plot=FALSE)
        hd<-hh$density
        mv1 <- 1:round(maxrv/2)
        mv2 <- (round(maxrv/2)+1):maxrv
        m1 <- mv1[which.max(hd[mv1 + 1])]
        m2 <- mv2[which.max(hd[mv2 + 1])]
        i1 <- which.min(hd[m1:m2 + 1]) - 1
        noiseFactor <- maxrv / i1
    }        
    rv[rv< maxrv/noiseFactor] <- 0
    cutVal <- maxrv/noiseFactor
    
    if(extreme)
        rv[rv > 0] <- 255

    ##rvals <- t(matrix(rv,nrow=nr,ncol=nc,byrow=TRUE))
    if(!onlyOne){
        ## reshape2::melt
        mrval <- cbind(expand.grid(1:nr,1:nc),Val=as.vector(rv))

        mx <- apply(rv,1,mean)
        #mx[mx<5] <- 0
        my <- apply(rv,2,mean)
        #my[my<5] <- 0
        difx <- which(diff(which(mx==0))>minPixelDiff)
        dify <- which(diff(which(my==0))>minPixelDiff)
        checkx <- length(difx)>1
        checky <- length(dify)>1

        if(any(c(checkx,checky))){
            nclust <- 2            
        }else{
            nclust <- 1
        }
        km <- stats::kmeans(mrval[mrval[,3]>0,1:2],nclust)
        
    }else{
        km <- list(cluster = rep(1,length(rv[rv>0])))
        nclust <- 1
    }
    km$cluster <- factor(km$cluster)
    levels(km$cluster) <- c("Left","Right")[order(km$centers[,1])]
    if(!assignSinglesByPosition & nclust == 1){
        levels(km$cluster) <- c("Unknown")
    }
    if(nlevels(km$cluster) > 1)
        km$cluster <- stats::relevel(km$cluster,"Left")

    res <- list()
    for(i in 1:nlevels(km$cluster)){
        rv3 <- rv
        rv3[rv3>0][as.integer(km$cluster) != i] <- 0
        if(pixelwise){
            cc <- Conte((rv3))
            cont <- cbind(cc$X,cc$Y)
        }else{
            cc<-contourLines(1:nr,1:nc,rv3,levels=ifelse(extreme,255,cutVal))
            cl <- unlist(lapply(cc,function(x)sum(sapply(2:length(x$x),
                                                      function(i)LineLength(c(x$x[i],x$y[i]),
                                                                            c(x$x[i-1],x$y[i-1]))))
                                ))
            ## cl <- unlist(lapply(cc,function(x)length(x$x)))
            cont <- do.call("cbind",cc[[which.max(cl)]][c("x","y")])
        }
        res[[i]] <- newStart_image(cont)
        attr(res[[i]],"Position") <- levels(km$cluster)[i]
        attr(res[[i]],"File") <- file
        attr(res[[i]],"NoiseFactor") <- noiseFactor
        attr(res[[i]],"ImagePixels") <- c(nc,nr)
        attr(res[[i]],"Normalized") <- FALSE
        attr(res[[i]],"Flipped") <- FALSE
        class(res[[i]]) <- "otolith_image"
    }
    return(res)
    
}

##' Normalize an Otolith Image
##'
##' @param dat 
##' @param n 
##' @param datCompare 
##' @param forceFlip 
##' @param flipByPosition 
##' @return A normalized Otolith Image
##' @author Christoffer Moesgaard Albertsen
##' @export
normalize_image <- function(dat,n,datCompare=NULL,forceFlip=FALSE,flipByPosition = c("No","Left","Right")){
    flipByPosition <- match.arg(flipByPosition)
    datNew <- rotate_image(dat)
    datNew <- center_image(datNew)
    datNew <- scale_image(datNew)
    datNew <- newStart_image(datNew)
    datNew <- setNumObs_image(datNew,n)
    if(!is.null(datCompare)){
        dcNorm <- normalize_image(datCompare,n)
        datNew <- flip_image(datNew,dcNorm,forceFlip)
    }
    if(flipByPosition == attr(dat,"Position"))
        datNew <- flip_image(datNew,datNew,TRUE)
    
    return(datNew)
}


