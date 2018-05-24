
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

flip_image <- function(dat,datCompare,forceFlip=FALSE, coord = 2){
    # Both dat and datCompare must be preprocessed
    dat1 <- dat
    dat2 <- dat
    mid <- mean(range(dat2[,coord]))
    dat2[,coord] <- mid - (dat2[,coord] - mid)
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
##' @importFrom utils tail
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
            imOut <- round((0.2989 * im[,,1] + 0.5870 * im[,,2] + 0.1140 * im[,,3]) * 255)
    }else{
        imOut <- round(im * 255)
    }
    return(imOut)
}


transformPixelMatrix <- function(pic,
                                 logisticTransform = FALSE,
                                 logisticTransformLocation = 255/2,
                                 logisticTransformScale = 1,
                                 gaussianBlur = FALSE,
                                 gaussianBlurSize = 3,
                                 gaussianBlurSigma = 1,
                                 unsharp = FALSE) {
    if(logisticTransform){
        pic <- plogis(pic,logisticTransformLocation,logisticTransformScale) * 255     
    }
    if(gaussianBlur & !unsharp){
            nr <- nrow(pic)
            nc <- ncol(pic)
            picOld <- pic
            borderCol <- mean(c(pic[c(1:(nr*0.02),nr:(nr-nr*0.02)),],pic[,c(1:(nc*0.02),nc:(nc-nc*0.02))]))
            pic <- matrix(borderCol,nr + (gaussianBlurSize-1),nc + (gaussianBlurSize-1))
            pic[((gaussianBlurSize-1)/2+1):(nrow(pic)-(gaussianBlurSize-1)/2),((gaussianBlurSize-1)/2+1):(ncol(pic)-(gaussianBlurSize-1)/2)] <- picOld
            pic <- gaussian_blur(pic, gaussianBlurSize, gaussianBlurSigma)
    }
    if(unsharp){
        nr <- nrow(pic)
        nc <- ncol(pic)
        picOld <- pic
        borderCol <- mean(c(pic[c(1:(nr*0.02),nr:(nr-nr*0.02)),],pic[,c(1:(nc*0.02),nc:(nc-nc*0.02))]))
        pic <- matrix(borderCol,nr + (gaussianBlurSize-1),nc+(gaussianBlurSize-1))
        pic[((gaussianBlurSize-1)/2+1):(nrow(pic)-(gaussianBlurSize-1)/2),((gaussianBlurSize-1)/2+1):(ncol(pic)-(gaussianBlurSize-1)/2)] <- picOld
        pic <- round(unsharpmask(pic, gaussianBlurSize, gaussianBlurSigma))
        pic[pic < 0] <- 0
        pic[pic > 255] <- 255
    }
    return(pic)
}

##' Read Otolith Images and Extract Contours
##'
##' @param file Image file path
##' @param noiseFactor Scalar value determining cutoff value for background noise on image. If NULL, a value is calculated automatically.
##' @param onlyOne Boolean value. Is there only one otolith on the image?
##' @param minPixelDiff Minimum pixel difference between otoliths
##' @param extreme Boolean value. Should pixel values be converted to 0/1?
##' @param pixelwise Boolean value. If TRUE, a pixel-wise algorithm is used to extract contours; otherwise, \code{grDevices::contourLines} is used.
##' @param assignSinglesByPosition Should single otoliths be assigned to Left/Right based on position on image?
##' @param minCountScale See details
##' @return otolith image information
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats kmeans relevel
##' @importFrom graphics hist
##' @importFrom grDevices contourLines
##' @export
read_image<- function(file,
                      noiseFactor = NULL,
                      onlyOne = FALSE,
                      minPixelDiff = 0.05 * min(nc,nr),
                      extreme = TRUE,
                      borderBasedCutOff = FALSE,
                      logisticTransform = FALSE,
                      logisticTransformLocation = c("mean","median","borderMean","borderMedian"),
                      logisticTransformScale = 1,
                      gaussianBlur = FALSE,
                      gaussianBlurSize = round(0.02 * min(nc,nr)),
                      gaussianBlurSigma = gaussianBlurSize / 10,
                      unsharp = FALSE,
                      pixelwise = FALSE,
                      assignSinglesByPosition = TRUE,
                      minCountScale = 0,
                      minCountForMax = 1e-4,
                      zeroCutOffPercent = 0.05){
    r<-getPixelMatrix(file)
    rv <- t(r[nrow(r):1,])
    nc <- ncol(rv)
    nr <- nrow(rv)

    
    whiteBorder <- mean(c(rv[c(1:(nr*0.02),nr:(nr-nr*0.02)),],rv[,c(1:(nc*0.02),nc:(nc-nc*0.02))])) > 255/2
    if(whiteBorder){
        rv <- 255 - rv
    }

    if(!is.numeric(logisticTransformLocation)){
        logisticTransformLocation <- match.arg(logisticTransformLocation)
    }
    ltl <- NA
    if(logisticTransform){
        bordervec <- c(rv[c(1:(nr*0.02),nr:(nr-nr*0.02)),],rv[,c(1:(nc*0.02),nc:(nc-nc*0.02))])
        if(is.numeric(logisticTransformLocation)){
            ltl <- logisticTransformLocation[1]
        }else{
            ltl <- switch(logisticTransformLocation,
                          mean = mean(rv),
                          median = median(rv),
                          borderMean = mean(bordervec),
                          borderMedian = median(bordervec)
                          )
        }       
    }
    rv <- transformPixelMatrix(rv,
                               logisticTransform = logisticTransform,
                               logisticTransformLocation = ltl,
                               logisticTransformScale = logisticTransformScale,
                               gaussianBlur = gaussianBlur,
                               gaussianBlurSize = gaussianBlurSize,
                               gaussianBlurSigma = gaussianBlurSigma,
                               unsharp = unsharp)
    
    maxrv <- max(rv)

    if(is.null(noiseFactor)){
        if(borderBasedCutOff){
            i1 <- quantile(c(rv[c(1:(nr*0.02),nr:(nr-nr*0.02)),],rv[,c(1:(nc*0.02),nc:(nc-nc*0.02))]),0.99)
        }else{
            hh<-graphics::hist(rv,breaks=(-1):maxrv + 0.5,plot=FALSE)
            hd<-hh$density * (hh$counts > nc*nr * minCountScale)
            ## mv1 <- 1:round(maxrv/2)
            ## mv2 <- (round(maxrv/2)+1):maxrv
            mx <- max(which(hh$counts > nc*nr * minCountForMax))
            mv1 <- 1:round(mx/2)
            mv2 <- (round(mx/2)+1):mx
            m1 <- mv1[which.max(hd[mv1 + 1])]
            m2 <- mv2[which.max(hd[mv2 + 1])]
            i1 <- max(m1) + which.min(hd[m1:m2 + 1]) - 1
        }
        noiseFactor <- maxrv / i1 
    }        
    rv[rv< maxrv/noiseFactor] <- 0
    cutVal <- maxrv/noiseFactor
    
    if(extreme)
        rv[rv > 0] <- 255
    
    ##rvals <- t(matrix(rv,nrow=nr,ncol=nc,byrow=TRUE))
    lvls <- c("Left","Right")
    orderIndx <- 1
    if(!onlyOne){
        ## reshape2::melt
        mrval <- cbind(expand.grid(1:nr,1:nc),Val=as.vector(rv))

        mx <- apply(rv,1,mean)
        ##mx[mx<5] <- 0
        my <- apply(rv,2,mean)
        ##my[my<5] <- 0
        difx <- which(diff(which(mx < maxrv*zeroCutOffPercent))>minPixelDiff)
        dify <- which(diff(which(my < maxrv*zeroCutOffPercent))>minPixelDiff)
        checkx <- length(difx)>1
        checky <- length(dify)>1

        if(any(c(checkx,checky))){
            nclust <- 2
            if(!checkx & checky)
                orderIndx <- 2
        }else{
            nclust <- 1
        }
        km <- stats::kmeans(mrval[mrval[,3]>0,1:2],nclust)
        
    }else{
        km <- list(cluster = rep(1,length(rv[rv>0])))
        nclust <- 1
    }
    km$cluster <- factor(km$cluster)
    levels(km$cluster) <- c("Left","Right")[order(km$centers[,orderIndx])]
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
            cc<-grDevices::contourLines(1:nr,1:nc,rv3,levels=ifelse(extreme,255,cutVal))
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
        attr(res[[i]],"LogisticTransform") <- list(used = logisticTransform,
                                                   location = ltl,
                                                   scale = logisticTransformScale)
        attr(res[[i]],"GaussianBlur") <- list(used = gaussianBlur & !unsharp,
                                              size = gaussianBlurSize,
                                              sigma = gaussianBlurSigma)
        attr(res[[i]],"UnsharpMask") <- list(used = unsharp)
        attr(res[[i]],"ImagePixels") <- c(nc,nr)
        attr(res[[i]],"Normalized") <- FALSE
        attr(res[[i]],"Flipped") <- FALSE
        class(res[[i]]) <- "otolith_contour"
    }
    attr(res,"File") <- file
    attr(res,"NoiseFactor") <- noiseFactor
    attr(res,"LogisticTransform") <- list(used = logisticTransform,
                                          location = ltl,
                                          scale = logisticTransformScale)
    attr(res,"GaussianBlur") <- list(used = gaussianBlur & !unsharp,
                                     size = gaussianBlurSize,
                                     sigma = gaussianBlurSigma)
    attr(res,"UnsharpMask") <- list(used = unsharp)
    attr(res,"ImagePixels") <- c(nc,nr)
    class(res) <- "otolith_image"
    return(res)
    
}

##' Normalize an Otolith Image
##'
##' @param dat otolith_image object
##' @param n Number of coordinates to output
##' @param datCompare Data to compare with
##' @param forceFlip Force a flip of the otolith?
##' @param flipByPosition Position (Left/Right) to flip. Use 'No' if otoliths should not be flipped by position.
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
    attr(datNew,"Normalized") <- TRUE
    return(datNew)
}


