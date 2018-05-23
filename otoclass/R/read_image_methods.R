
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
    cat("Position:",attr(x,"Position"))
    invisible(x)
}

##' Print otolith image
##'
##' @param x Otolith image
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method print otolith_contour
##' @export
print.otolith_image <- function(x, ...){
    cat("Otolith image of", attr(x,"File"),"\n")
    cat("Containing",length(x),ifelse(length(x)==1,"contour.","contours."))
    invisible(x)
}


##' Plot otolith image
##'
##' @param x Otolith image
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method plot otolith_image
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
                                    unsharp = attr(x,"UnsharpMask")$size
                                    )
    }
    par(mar = c(0,0,0,0), oma = c(0,0,0,0))
    plot(0,0,xlim=c(1,ncol(pic)),ylim=c(1,nrow(pic)),asp=1, type = "n")
    rasterImage(pic / 255,1,1,ncol(pic),nrow(pic))
    invisible(lapply(x,lines, ...))
    invisible(x)
}


##' Plot otolith contour
##'
##' @param x Otolith contour
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method plot otolith_contour
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
                                    unsharp = attr(x,"UnsharpMask")$used
                                    )
    }
    par(mar = c(0,0,0,0), oma = c(0,0,0,0))
    plot(0,0,xlim=c(1,ncol(pic)),ylim=c(1,nrow(pic)),asp=1, type = "n")
    rasterImage(pic / 255,1,1,ncol(pic),nrow(pic))
    lines(x, ...)
    invisible(x)
}


##' Plot otolith contour list
##'
##' @param x Otolith contour list
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method plot otolith_contour_list
##' @export
plot.otolith_contour_list<- function(x,asUsed = FALSE, ask = TRUE, ...){
    oldAsk <- par("ask")
    par(ask = ask)
    for(i in 1:length(x))
        plot(x[[i]], asUsed = asUsed, ...)
    par(ask = oldAsk)
}

