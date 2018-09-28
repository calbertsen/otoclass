
##' Rotate an otolith
##'
##' @param x object to rotate
##' @param angle Angles to rotate the object
##' @param degrees Is the angle in degrees?
##' @param ... not used
##' @return A rotated element 
##' @author Christoffer Moesgaard Albertsen
##' @export
rotate <- function(x, angle, degrees = TRUE, ...){
    UseMethod("rotate")
}

##' @export
rotate.matrix <- function(x, angle, degrees = TRUE, center = FALSE, newStart = FALSE, ...){
    if(center){
        cent <- apply(x,2,function(xx) mean(range(xx)))
        xc <- t(t(x) - cent)
    }else{
        xc <- x
    }
    if(degrees)
        angle <- angle / 180 * pi
    rotMat <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), ncol=2)
    val <- t( t(xc %*% rotMat) + cent)
    if(newStart)
        val <- newStart_image(val)
    return(val)
}

##' @export
rotate.otolith_contour <- function(x, angle, degrees = TRUE, ...){
    at <- attributes(x)
    val <- rotate(unclass(x), angle = angle, degrees = degrees, center = TRUE, newStart = TRUE, ...)
    attributes(val) <- at
    if(degrees)
        angle <- angle / 180 * pi
    attr(val, "Rotated") <- (attr(val, "Rotated") + angle) %% (2 * pi)
    return(val)
}


##' @export
rotate.otolith_contour_list <- function(x, angle, degrees = TRUE, ...){
    at <- attributes(x)
    val <- lapply(x, rotate, angle = angle, degrees = degrees, ...)
    attributes(val) <- at
    return(val)
}

##' @export
rotate.otolith_contour_image <- function(x, angle, degrees = TRUE, ...){
    at <- attributes(x)
    val <- lapply(x, rotate, angle = angle, degrees = degrees, ...)
    attributes(val) <- at
    return(val)
}
