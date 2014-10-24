
setNumObs_image<- function(dat,t){
    tmpX <- approx(x=seq(0,1,len=dim(dat)[1]),y=dat[,1],xout=t)
    tmpY <- approx(x=seq(0,1,len=dim(dat)[1]),y=dat[,2],xout=t)
 
    return(cbind(tmpX$y,tmpY$y))
}
center_image <- function(dat){
    tmpX <- dat[,1]
    tmpY <- dat[,2]
    return(cbind(tmpX-mean(tmpX),
                 tmpY-mean(tmpY))
           )
}
scale_image <- function(dat){
    tmpX <- dat[,1]
    tmpY <- dat[,2]
    return(cbind(tmpX/diff(range(tmpX)),
                 tmpY/diff(range(tmpX)))
           )

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
    return(t(apply(val,1,function(x)x+oldCenter)))
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
    return(datNew)
}
flip_image <- function(dat,datCompare,forceFlip=FALSE){
    # Both dat and datCompare must be preprocessed
    dat1 <- dat
    dat2 <- dat
    dat2[,2] <- -dat2[,2]
    dat2 <-  newStart_image(dat2)
    mse1 <- mean((dat1-datCompare)^2)
    mse2 <- mean((dat2-datCompare)^2)

    if(forceFlip)
        return(dat2)
    
    if(mse1 <= mse2){
        return(dat1)
    }else{
        return(dat2)
    }
}

#' @export
read_image<- function(file,noiseFactor = 4, onlyOne = FALSE, minPixelDiff = 20){
    require(raster)
    library(reshape2)
    r<-raster(file)
    r[r[]< max(r[])/noiseFactor] <- 0

    if(!onlyOne){
        rvals <- t(as.matrix(r))
        mrval <- melt(rvals)

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
        km <- kmeans(mrval[mrval[,3]>0,1:2],nclust)

                                        #The real time problem
    }else{
        km <- list(cluster = rep(1,length(r[r>0])))
        nclust <- 1
    }

    res <- list()
    for(i in 1:nclust){
        r3 <- r
        r3[r3>0][km$cluster!=i] <- 0 # 1 sec
        c<-rasterToContour(r3) # 0.35 sec
        cl <- unlist(lapply(c@lines[[1]]@Lines,function(x)dim(x@coords)[1]))
        res[[i]] <- c@lines[[1]]@Lines[[which.max(cl)]]@coords
    }

    return(res)
    
}

#' @export
normalize_image <- function(dat,t,datCompare=NULL,forceFlip=FALSE){
    datNew <- rotate_image(dat)
    datNew <- center_image(datNew)
    datNew <- scale_image(datNew)
    datNew <- newStart_image(datNew)
    if(!is.null(datCompare)){
        datNew <- flip_image(datNew,datCompare,forceFlip)
    }   
    datNew <- setNumObs_image(datNew,t)
    return(datNew)
}


