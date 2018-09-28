##' Interpolate otolith coordinates through EFDs
##'
##' @param x otolith object
##' @param N Number of coordinates to output
##' @param Nefd Number of EFDs to use
##' @param ... not used
##' @return an interpolated otolith object
##' @author Christoffer Moesgaard Albertsen
##' @export
interpolate <- function(x, N, Nefd, ...){
    UseMethod("interpolate")
}

##' @export
interpolate.matrix <- function(x, N, Nefd, ...){
    eList <- efd(x, N=N, returnAsList = TRUE)
    eList$n <- N
    val <- do.call(efd2coord, eList)
    #val <-  newStart_image(val)
    return(val)
}

##' @export
interpolate.otolith_contour <- function(x, N, Nefd, ...){
    at <- attributes(x)
    val <- interpolate(unclass(x), N=N, Nefd = Nefd, ...)
    at$dim[1] <- N
    attributes(val) <- at
    return(val)
}


##' @export
interpolate.otolith_contour_list <- function(x, N, Nefd, ...){
    at <- attributes(x)
    val <- lapply(x, interpolate, N=N, Nefd = Nefd, ...)
    attributes(val) <- at
    return(val)
}

##' @export
interpolate.otolith_contour_image<- function(x, N, Nefd, ...){
    at <- attributes(x)
    val <- lapply(x, interpolate, N=N, Nefd = Nefd, ...)
    attributes(val) <- at
    return(val)
}
