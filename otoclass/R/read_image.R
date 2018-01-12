
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

##' Read Otolith Images and Extract Contours
##'
##' @param file 
##' @param noiseFactor 
##' @param onlyOne 
##' @param minPixelDiff 
##' @return otolith image information
##' @author Christoffer Moesgaard Albertsen
##' @importFrom raster raster as.matrix rasterToContour ncol nrow
##' @importFrom stats kmeans relevel
##' @importFrom reshape2 melt
##' @export
read_image<- function(file,noiseFactor = NULL, onlyOne = FALSE, minPixelDiff = 0.03 * min(nc,nr), extreme = FALSE, assignSinglesByPosition = TRUE){
    r<-raster::raster(file)
    nc <- raster::ncol(r)
    nr <- raster::nrow(r)

    whiteBorder <- mean(c(r[c(1:(nr*0.02),nr:(nr-nr*0.02)),],r[,c(1:(nc*0.02),nc:(nc-nc*0.02))])) > 255/2
    if(whiteBorder){
        r[] <- 255 - r[]
    }

    if(is.null(noiseFactor)){
        hh<-hist(r[],breaks=0:255,plot=FALSE)
        hd<-hh$density
        i1 <- which.min(hd[1:length(hd) > which.max(hd)]) + which.max(hd) - 1
        noiseFactor <- max(r[]) / i1
        }        
    r[r[]< max(r[])/noiseFactor] <- 0
    
    if(extreme)
        r[r[] > 0] <- 255

    if(!onlyOne){
        rvals <- t(raster::as.matrix(r))
        mrval <- reshape2::melt(rvals)

        mx <- apply(rvals,1,mean)
        #mx[mx<5] <- 0
        my <- apply(rvals,2,mean)
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

                                        #Can perhaps be faster? 0.4 sec
        km <- stats::kmeans(mrval[mrval[,3]>0,1:2],nclust)
        

                                        #The real time problem
    }else{
        km <- list(cluster = rep(1,length(r[r>0])))
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
        r3 <- r
        r3[r3>0][as.integer(km$cluster) != i] <- 0 # 1 sec
        if(extreme){
            cc <- Conte(raster::as.matrix(r3))
            cont <- cbind(cc$X,cc$Y)
        }else{
            cc<-raster::rasterToContour(r3) # 0.35 sec
            cl <- unlist(lapply(cc@lines[[1]]@Lines,function(x)dim(x@coords)[1]))
            cont <- cc@lines[[1]]@Lines[[which.max(cl)]]@coords
        }
        res[[i]] <- cont
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


