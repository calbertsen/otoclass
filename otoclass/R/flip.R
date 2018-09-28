##' Flip an otolith
##'
##' @param x otolith object
##' @param axis axis to flip around
##' @param ... not used
##' @return a flipped otolith object
##' @author Christoffer Moesgaard Albertsen
##' @export
flip <- function(x, axis, ...){
    UseMethod("flip")
}

##' @export
flip.matrix <- function(x, axis, ...){
    if(length(axis) != 1 || !match(axis,1:2,FALSE))
        stop("axis must be 1 or 2")
    val <- x
    cent <- mean(range(val[,axis]))
    val[,axis] <- cent - (val[,axis] - cent)
    val <-  newStart_image(val)
    return(val)
}

##' @export
flip.otolith_contour <- function(x, axis, ...){
    at <- attributes(x)
    val <- flip(unclass(x), axis = axis, ...)
    attributes(val) <- at
    attr(val,"Flipped")[axis] <- !attr(val,"Flipped")[axis]
    return(val)
}


##' @export
flip.otolith_contour_list <- function(x, axis, ...){
    at <- attributes(x)
    val <- lapply(x, flip, axis = axis, ...)
    attributes(val) <- at
    return(val)
}

##' @export
flip.otolith_contour_image<- function(x, axis, ...){
    at <- attributes(x)
    val <- lapply(x, flip, axis = axis, ...)
    attributes(val) <- at
    return(val)
}
