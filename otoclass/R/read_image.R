
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

newStart_image <- function(dat, polar = FALSE){
    isClosed <- isTRUE(all.equal(dat[1,],dat[nrow(dat),]))
    if(isClosed)
        dat <- dat[-nrow(dat),]
    ldat <- nrow(dat)
    if(polar){
        newBegin <- which.min(Arg(complex(re=dat[,1],im=dat[,2])))
    }else{
        newBegin <- which.max(dat[,1])
    }
    if(newBegin == 1){
        datNew <- dat
    }else{
        datNew <- rbind(dat[newBegin:ldat,],
                        dat[1:(newBegin-1),]
                        )
    }
    if(datNew[2,2] < datNew[1,2])
        datNew <- apply(datNew,2,rev)
    if(isClosed){
        datNew <- rbind(datNew,datNew[1,])
    }
    attrs <- attributes(dat)
    aN <- attributes(datNew)
    attrs[intersect(names(aN),names(attrs))] <- attributes(datNew)[intersect(names(aN),names(attrs))]
    attributes(datNew) <- attrs
    return(datNew)
}

flip_image <- function(dat,datCompare,forceFlip=FALSE, coord = 2){
                                        # Both dat and datCompare must be preprocessed
    attri <- attributes(dat)
    dat1 <- dat
    dat2 <- dat
    mid <- mean(range(dat2[,coord]))
    dat2[,coord] <- mid - (dat2[,coord] - mid)
    dat2 <-  newStart_image(dat2)
    attri$Flipped[coord] <- attri$Flipped[coord]
    attributes(dat2) <- attri
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

extendPic <- function(pic, kernelSize){
    nr <- nrow(pic)
    nc <- ncol(pic)
    picOld <- pic
    borderCol <- mean(c(pic[c(1:(nr*0.02),nr:(nr-nr*0.02)),],pic[,c(1:(nc*0.02),nc:(nc-nc*0.02))]))
    pic <- matrix(borderCol,nr + (kernelSize-1),nc + (kernelSize-1))
    pic[((kernelSize-1)/2+1):(nrow(pic)-(kernelSize-1)/2),((kernelSize-1)/2+1):(ncol(pic)-(kernelSize-1)/2)] <- picOld
    pic
}

whiteBalance <- function(pic, pct = 0.00025){
    rng <- quantile(pic,c(pct,1-pct))
    newPic <- (pic - rng[1]) * 255 / diff(rng)
    newPic[newPic < 0] <- 0
    newPic[newPic > 255] <- 255
    newPic
}

bruteCutCol <- function(pic, pct = 0.05){
    cm <- apply(pic,2,max) #colMeans(pic)
    rm <- apply(pic,1,max) #rowMeans(pic)
    cutValC <- quantile(cm[seq(floor(ncol(pic)/2 -  ncol(pic)*pct),
                               ceiling(ncol(pic)/2 +  ncol(pic)*pct), 1)],0.75)
     ## Cut columns
    pic[,cm < cutValC] <- 0
    ii <- which(cm >= cutValC)
    iChunks <- split(ii, cumsum(c(0,diff(ii)>1)))
    iCL <- sapply(iChunks,length)
    pic[,unlist(iChunks[iCL < min(tail(sort(iCL),2))])] <- 0
    ## Cut rows
    pic[rm < cutValC,] <- 0
    ii <- which(rm >= cutValC)
    iChunks <- split(ii, cumsum(c(0,diff(ii)>1)))
    iCL <- sapply(iChunks,length)
    pic[unlist(iChunks[iCL < min(tail(sort(iCL),1))]),] <- 0
    #pic[pic < 0.5 * cutValC] <- 0
    pic
}

histogramEqualization <- function(pic){
    pn <- table(factor(pic,0:255)) / length(pic)
    Pn <- cumsum(pn)
    newPic <- matrix(floor(255 * Pn[pic+1]), nrow(pic), ncol(pic))
    newPic
}

blending <- function(pic, alpha, FUN = histogramEqualization, ...){
    FUN(pic,...) * alpha + (1 - alpha) * pic
}

histogramCut <- function(pic, alpha){
    pic[pic < quantile(pic,alpha)] <- 0
    pic
}

contrastLimitedHistogramEqualization <- function(pic, C = 10){
    pn0 <- table(factor(pic,0:255)) / length(pic)
    ## Find actual clipping limit
    C <- 100
    pn0N <- pn0 * length(pic)
    CN <- C * mean(pn0N) 
    top <- CN 
    bottom <- 0
    while(top - bottom > 1){
        middle <- (top + bottom) / 2
        excess <- (pn0N - middle) * as.numeric(pn0N > middle)
        S <- sum(excess)
        if(S > (CN-middle) * 256){
            top <- middle
        }else{
            bottom <- middle
        }
    }
    ## Clip
    pn <- pmin(pn0N + 1, bottom)
    ## Re-normalize
    pn <- pn / sum(pn)
    ## Equalize
    Pn <- cumsum(pn)
    newPic <- matrix(floor(255 * Pn[pic+1]), nrow(pic), ncol(pic))
    newPic
}

##' @importFrom stats plogis
transformPixelMatrix <- function(pic,
                                 transformations = character(0),
                                 logisticTransformLocation = 255/2,
                                 logisticTransformScale = 1,
                                 gaussianBlurSize = 3,
                                 gaussianBlurSigma = 1,
                                 floodFillTolerance = 0.1,
                                 floodFillCol = 0,
                                 whiteBalance = 0.25,
                                 alpha = 0.4) {
    if(length(transformations) == 0)
        return(pic)
    tIndx <- pmatch(transformations, c("logistic","gaussianBlur","unsharp","floodFill","sobel","edgeDetect","centerFloodFill","whiteBalance","histogramEqualization","histogramCut","bruteCutCol","medianfilter","meanfilter"), duplicates.ok = TRUE)
    ## Size of input
    resizeInput <- function(y,i){
        xx <- rep(NA, length(tIndx))
        xx[tIndx %in% i] <- rep(y, length.out = sum(tIndx %in% i))
        xx
    }
    logisticTransformLocation = resizeInput(logisticTransformLocation,1)
    logisticTransformScale = resizeInput(logisticTransformScale,1)
    gaussianBlurSize = resizeInput(gaussianBlurSize, c(2,3,12,13))
    gaussianBlurSigma = resizeInput(gaussianBlurSigma, c(2,3))
    floodFillTolerance = resizeInput(floodFillTolerance, c(4,7))
    whiteBalance = resizeInput(whiteBalance, c(8))
    alpha = resizeInput(alpha, c(10,11))
    floodFillCol = resizeInput(floodFillCol, c(4,7))
    
    for(ii in seq_along(tIndx)){
        vv <- tIndx[ii]
        if(is.na(vv)){
            warning(sprintf("Wrong transformation option: %s", transformations[ii]))
        }else if(vv == 1){ ## logistic transform
            pic <- stats::plogis(pic,logisticTransformLocation[ii],
                                 logisticTransformScale[ii]) * 255     
        }else if(vv == 2){ ## Gaussian blur
               pic <- gaussian_blur(extendPic(pic, gaussianBlurSize[ii]), gaussianBlurSize[ii], gaussianBlurSigma[ii])
        }else if(vv == 3){ ## Unsharp          
            pic <- round(unsharpmask(extendPic(pic, gaussianBlurSize[ii]), gaussianBlurSize[ii], gaussianBlurSigma[ii]))
            pic[pic < 0] <- 0
            pic[pic > 255] <- 255
        }else if(vv == 4){ ## Flood fill
            pic <- floodfill(pic, floodFillTolerance[ii], 0, 0, TRUE, floodFillCol[ii])
        }else if(vv == 5){ ## Sobel           
            pic <- to01(sobel(extendPic(pic,3))) * 255
        }else if(vv == 6){ ## edgeDetect
            pic <- to01(edge_detect(extendPic(pic,3))) * 255
        }else if(vv == 7){ ## Center flood fill
            pic <- floodfill(to01(pic)*255, floodFillTolerance[ii], ncol(pic)/2, nrow(pic)/2, TRUE, floodFillCol[ii])
        }else if(vv == 8){ ## white balance
            pic <- whiteBalance(to01(pic)*255, whiteBalance[ii])
        }else if(vv == 9){ ## histogram equalization
            pic <- histogramEqualization(to01(pic)*255)
        }else if(vv == 10){ ## histogram cut
            pic <- histogramCut(to01(pic)*255, alpha[ii])
        }else if(vv == 11){ ## Brute cut columns
            pic <- t(bruteCutCol(to01(t(pic))*255, alpha[ii]))
        }else if(vv == 12){
            pic <- medianfilter(pic, gaussianBlurSize[ii])
        }else if(vv == 13){
            pic <- meanfilter(extendPic(pic,gaussianBlurSize[ii]), gaussianBlurSize[ii])
        }
    }
    return(round(pic))
}

##' Read Otolith Images and Extract Contours
##'
##' @param file Image file path
##' @param noiseFactor Scalar value determining cutoff value for background noise on image. If NULL, a value is calculated automatically.
##' @param onlyOne Boolean value. Is there only one otolith on the image?
##' @param minPixelDiff Minimum pixel difference between otoliths
##' @param extreme Boolean value. Should pixel values be converted to 0/1?
##' @param borderBasedCutOff Use a border based cut off?
##' @param transformations Vector of image transformations: logistic, gaussianBlur, unsharp, floodFill
##' @param logisticTransformLocation Location parameter for logistic transformation
##' @param logisticTransformScale Scale parameter for logistic transformation
##' @param gaussianBlurSize Size (n x n) for Gaussian blur kernel
##' @param gaussianBlurSigma Variance parameter of Gaussian distribution used to blur
##' @param pixelwise Boolean value. If TRUE, a pixel-wise algorithm is used to extract contours; otherwise, \code{grDevices::contourLines} is used.
##' @param assignSinglesByPosition Should single otoliths be assigned to Left/Right based on position on image?
##' @param minCountScale See details
##' @param minCountForMax See details
##' @param reduceCutOffPercent See details
##' @param zeroCutOffPercent See details
##' @return otolith image object
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats kmeans relevel median quantile
##' @importFrom graphics hist
##' @importFrom grDevices contourLines
##' @export
read_image<- function(file,
                      noiseFactor = NULL,
                      onlyOne = FALSE,
                      minPixelDiff = 0.05 * min(nc,nr),
                      extreme = TRUE,
                      floodFillTolerance = 0.1,
                      floodFillCol = 0,
                      whiteBalance = 0.4,
                      histogramAlpha = 0.4,
                      borderBasedCutOff = FALSE,
                      transformations = character(0),
                      logisticTransformLocation = c("mean","median","borderMean","borderMedian"),
                      logisticTransformScale = 1,
                      gaussianBlurSize = round(0.02 * min(nc,nr)),
                      gaussianBlurSigma = gaussianBlurSize / 10,
                      pixelwise = FALSE,
                      assignSinglesByPosition = TRUE,
                      forceCenter = FALSE,
                      minCountScale = 0,
                      minCountForMax = 1e-4,
                      reduceCutOffPercent = 0,
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
        bordervec <- c(rv[c(1:(nr*0.02),nr:(nr-nr*0.02)),],rv[,c(1:(nc*0.02),nc:(nc-nc*0.02))])
        logisticTransformLocation <- sapply(logisticTransformLocation, function(ltf){
            switch(ltf,
                   mean = mean(rv),
                   median = stats::median(rv),
                   borderMean = mean(bordervec),
                   borderMedian = stats::median(bordervec)
                   )
        })       
    }
    rv <- transformPixelMatrix(rv,
                               transformations = transformations,
                               logisticTransformLocation = logisticTransformLocation,
                               logisticTransformScale = logisticTransformScale,
                               gaussianBlurSize = gaussianBlurSize,
                               gaussianBlurSigma = gaussianBlurSigma,
                               floodFillTolerance = floodFillTolerance,
                               whiteBalance = whiteBalance,
                               alpha = histogramAlpha,
                               floodFillCol = floodFillCol)
   
    maxrv <- max(rv)

    if(is.null(noiseFactor)){
        if(borderBasedCutOff){
            i1 <- stats::quantile(c(rv[c(1:(nr*0.02),nr:(nr-nr*0.02)),],rv[,c(1:(nc*0.02),nc:(nc-nc*0.02))]),0.99)
        }else{
            hh<-graphics::hist(rv,breaks=(-1):maxrv + 0.5,plot=FALSE)
            hd<-hh$density * (hh$counts > nc*nr * minCountScale)
            ## mv1 <- 1:round(maxrv/2)
            ## mv2 <- (round(maxrv/2)+1):maxrv
            mx <- max(which(hh$counts > nc*nr * minCountForMax))
            mv1 <- 2:round(mx/2)
            mv2 <- (round(mx/2)+1):mx
            m1 <- mv1[which.max(hd[mv1 + 1])]
            m2 <- mv2[which.max(hd[mv2 + 1])]
            i1 <- max(m1) + which.min(hd[m1:m2 + 1]) - 1
            i1 <- i1 - (i1 - m1) * reduceCutOffPercent
        }
        noiseFactor <- maxrv / i1 
    }
    if(!is.na(noiseFactor)){
        rv[rv< maxrv/noiseFactor] <- 0
        cutVal <- maxrv/noiseFactor
    }else{
        cutVal <- NA
    }
    
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
        km <- list(cluster = factor(rep(1,length(rv[rv>0]))),
                   centers = matrix(c(nr/2,nc/2),2,1))
        nclust <- 1
    }
    km$cluster <- factor(km$cluster)
    if(onlyOne){
        levels(km$cluster) <- c("Unknown")
    }else{
        levels(km$cluster) <- c("Left","Right")[order(km$centers[,orderIndx])]
    }
    if(assignSinglesByPosition & nclust == 1){
        if(km$centers[1,1] < nr/2){
            levels(km$cluster) <- c("Left")
        }else{
            levels(km$cluster) <- c("Right")
        }
    }
    if(!assignSinglesByPosition & nclust == 1){
        levels(km$cluster) <- c("Unknown")
    }
    if(nlevels(km$cluster) > 1)
        km$cluster <- stats::relevel(km$cluster,"Left")

    res <- list()
    for(i in 1:nlevels(km$cluster)){
        rv3 <- rv
        rv3[rv3>0][as.integer(km$cluster) != i] <- NA
        if(pixelwise){
            cc <- Conte((rv3))
            cont <- cbind(cc$X,cc$Y)
        }else{
            if(is.na(cutVal)){
                if(extreme){
                    cc<-grDevices::contourLines(1:nr,1:nc,rv3, levels = 255)
                }else{
                    cc<-grDevices::contourLines(1:nr,1:nc,rv3)
                }
            }else{
                cc<-grDevices::contourLines(1:nr,1:nc,rv3,levels=ifelse(extreme,255,cutVal))
            }
            if(forceCenter){
                hasCenter <- sapply(cc, function(X) point_in_polygon(px=nr/2,py=nc/2,X=cbind(X$x,X$y)))
            }else{
                hasCenter <- rep(TRUE, length(cc))
            }
            cc <- cc[hasCenter]
            cl <- unlist(lapply(cc,function(x)sum(sapply(2:length(x$x),
                                                         function(i)LineLength(c(x$x[i],x$y[i]),
                                                                               c(x$x[i-1],x$y[i-1]))))
                                ))
            ca <- sapply(cc,function(x)abs(otoclass:::polygon_area(cbind(x$x,x$y))))
            ## cl <- unlist(lapply(cc,function(x)length(x$x)))
            cont <- do.call("cbind",cc[[which.max(ca)]][c("x","y")])
        }
        res[[i]] <- newStart_image(cont)
        attr(res[[i]],"Position") <- levels(km$cluster)[i]
        attr(res[[i]],"File") <- file
        attr(res[[i]],"NoiseFactor") <- noiseFactor
        attr(res[[i]],"Transformation") <- list(transformations = transformations,
                                                logisticTransformLocation = logisticTransformLocation,
                                                logisticTransformScale = logisticTransformScale,
                                                gaussianBlurSize = gaussianBlurSize,
                                                gaussianBlurSigma = gaussianBlurSigma,
                                                floodFillTolerance = floodFillTolerance,
                                                floodFillCol = floodFillCol,
                                                whiteBalance = whiteBalance,
                                                alpha = histogramAlpha)
        attr(res[[i]],"ImagePixels") <- c(nc,nr)
        attr(res[[i]],"Normalized") <- FALSE
        attr(res[[i]],"Flipped") <- c(FALSE,FALSE)
        attr(res[[i]],"Rotated") <- 0
        attr(res[[i]],"Area") <- polygon_area(res[[i]])
        class(res[[i]]) <- "otolith_contour"
    }
    attr(res,"File") <- file
    attr(res,"NoiseFactor") <- noiseFactor
    attr(res,"Transformation") <- list(transformations = transformations,
                                       logisticTransformLocation = logisticTransformLocation,
                                       logisticTransformScale = logisticTransformScale,
                                       gaussianBlurSize = gaussianBlurSize,
                                       gaussianBlurSigma = gaussianBlurSigma,
                                       floodFillTolerance = floodFillTolerance,
                                       floodFillCol = floodFillCol,
                                       whiteBalance = whiteBalance,
                                       alpha = histogramAlpha)
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
normalize_image <- function(dat,n,datCompare=NULL,forceFlip=FALSE,flipByPosition = c("No","Left","Right", startByPolar = FALSE)){
    flipByPosition <- match.arg(flipByPosition)
    datNew <- rotate_image(dat)
    datNew <- center_image(datNew)
    datNew <- scale_image(datNew)
    datNew <- newStart_image(datNew, startByPolar)
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


