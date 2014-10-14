
#' @export

read_image <- function(file){
    require(raster)
    
    r<-raster(file)
    c<-rasterToContour(r)

    cl <- unlist(lapply(c@lines[[1]]@Lines,function(x)dim(x@coords)[1]))
    
    c@lines[[1]]@Lines[[which.max(cl)]]@coords
}

#' @export

normalize_image<- function(dat,t=seq(0,1,len=1000)){

    tmpX <- approx(x=seq(0,1,len=dim(dat)[1]),y=dat[,1],xout=t)
    tmpY <- approx(x=seq(0,1,len=dim(dat)[1]),y=dat[,2],xout=t)
 
   return(cbind((tmpX$y - mean(range(tmpX$y)))/abs(diff(range(tmpX$y))),
                     (tmpY$y - mean(range(tmpY$y)))/abs(diff(range(tmpX$y))))
          )

}
