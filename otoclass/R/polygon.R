##' Signed area of a closed polygon
##'
##' The function uses the shoelace formula to calculate the signed area of a closed non-intersecting polygon.
##' @param X A numeric matrix with two columns of vertices
##' @return Polygon area in same units as input
##' @author Christoffer Moesgaard Albertsen
polygon_area <- function(X){
    ## Check dims
    if(!is.numeric(X) | !is.matrix(X) | ncol(X) != 2 )
        stop("X must be a numeric matrix with two columns.")
    ## Check if closed contour
    if(!all(X[1,] == X[nrow(X),]))
        X <- rbind(X,X[1,])
    ## Calculate
    .Call("polygon_area", x=X[,1], y=X[,2], PACKAGE = "otoclass")
}

##' Centroid of a closed polygon
##'
##' The function calculates the centroid of a closed polygon.
##' @param X A numeric matrix with two columns of vertices
##' @return Polygon centroid in same units as input
##' @author Christoffer Moesgaard Albertsen
polygon_centroid <- function(X){
    ## Check dims
    if(!is.numeric(X) | !is.matrix(X) | ncol(X) != 2 )
        stop("X must be a numeric matrix with two columns.")
    ## Check if closed contour
    if(!all(X[1,] == X[nrow(X),]))
        X <- rbind(X,X[1,])
    ## Calculate
    .Call("polygon_centroid", x=X[,1], y=X[,2], PACKAGE = "otoclass")
}
