
##' Print otolith image contour
##'
##' @param x Otolith image contour
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method print otolith_contour
##' @export
print.otolith_contour <- function(x, ...){
    cat("Otolith contour of", attr(x,"File"),"\n")
    cat("Position:",attr(x,"Position"),"\n")
    invisible(x)
}

##' Print otolith image
##'
##' @param x Otolith image
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method print otolith_image
##' @export
print.otolith_image <- function(x, ...){
    cat("Otolith image of", attr(x,"File"),"\n")
    cat("Containing",length(x),ifelse(length(x)==1,"contour.","contours."),"\n")
    invisible(x)
}


##' Plot otolith image
##'
##' @param x Otolith image
##' @param asUsed Plot image as read; i.e. transformed
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method plot otolith_image
##' @importFrom graphics par plot rasterImage lines
##' @export
plot.otolith_image <- function(x,asUsed = FALSE, ...){
    pic <- getPixelMatrix(attr(x,"File"))
    if(asUsed){
        pic <- transformPixelMatrix(pic,
                                    logisticTransform = attr(x,"LogisticTransform")$used,
                                    logisticTransformLocation = attr(x,"LogisticTransform")$location,
                                    logisticTransformScale = attr(x,"LogisticTransform")$scale,
                                    gaussianBlur = attr(x,"GaussianBlur")$used,
                                    gaussianBlurSize = attr(x,"GaussianBlur")$size,
                                    gaussianBlurSigma = attr(x,"GaussianBlur")$sigma,
                                    unsharp = attr(x,"UnsharpMask")$used,
                                    floodFillTolerance = attr(x,"floodFillTolerance")
                                    )
    }
    graphics::par(mar = c(0,0,0,0), oma = c(0,0,0,0))
    graphics::plot(0,0,xlim=c(1,ncol(pic)),ylim=c(1,nrow(pic)),asp=1, type = "n")
    graphics::rasterImage(pic / 255,1,1,ncol(pic),nrow(pic))
    invisible(lapply(x,graphics::lines, ...))
    invisible(x)
}


##' Plot otolith contour
##'
##' @param x Otolith contour
##' @param asUsed Plot image as read; i.e. transformed
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method plot otolith_contour
##' @importFrom graphics par plot rasterImage lines
##' @export
plot.otolith_contour <- function(x,asUsed = FALSE, ...){
    pic <- getPixelMatrix(attr(x,"File"))
    if(asUsed){
        pic <- transformPixelMatrix(pic,
                                    logisticTransform = attr(x,"LogisticTransform")$used,
                                    logisticTransformLocation = attr(x,"LogisticTransform")$location,
                                    logisticTransformScale = attr(x,"LogisticTransform")$scale,
                                    gaussianBlur = attr(x,"GaussianBlur")$used,
                                    gaussianBlurSize = attr(x,"GaussianBlur")$size,
                                    gaussianBlurSigma = attr(x,"GaussianBlur")$sigma,
                                    unsharp = attr(x,"UnsharpMask")$used,
                                    floodFillTolerance = attr(x,"floodFillTolerance")
                                    )
    }
    graphics::par(mar = c(0,0,0,0), oma = c(0,0,0,0))
    graphics::plot(0,0,xlim=c(1,ncol(pic)),ylim=c(1,nrow(pic)),asp=1, type = "n")
    graphics::rasterImage(pic / 255,1,1,ncol(pic),nrow(pic))
    graphics::lines(x, ...)
    invisible(x)
}


##' Plot otolith contour list
##'
##' @param x Otolith contour list
##' @param asUsed Plot image as read; i.e. transformed
##' @param ask Ask before plotting next contour?
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method plot otolith_contour_list
##' @export
plot.otolith_contour_list<- function(x,asUsed = FALSE, ask = TRUE, ...){
    oldAsk <- graphics::par("ask")
    graphics::plot(x[[1]], asUsed = asUsed, ...)
    graphics::par(ask = ask)
    if(length(x) > 1)
        for(i in 2:length(x))
            graphics::plot(x[[i]], asUsed = asUsed, ...)
    graphics::par(ask = oldAsk)
}


##' @export
##' @author Christoffer Moesgaard Albertsen
`[.otolith_contour_list` <- function(x, i){
    a <- attributes(x)
    r <- unclass(x)[i]
    attributes(r) <- a
    return(r)
}

##' @export
##' @author Christoffer Moesgaard Albertsen
`[[<-.otolith_contour_list` <- function(x, i, value){
    a <- attributes(x)
    r <- unclass(x) ## Convert to list
    r[[i]] <- value
    attributes(r) <- a
    return(r)
}

##' @export
##' @author Christoffer Moesgaard Albertsen
`[<-.otolith_contour_list` <- function(x, i, value){
    a <- attributes(x)
    r <- unclass(x) ## Convert to list
    r[i] <- value
    attributes(r) <- a
    return(r)
}

