
##' Convert to otolith_contour_list
##' 
##' @param x otolith_contour object or output from read_image
##' @param ... other parameters
##' @return otolith_contour_list object
##' @author Christoffer Moesgaard Albertsen
as.otolith_contour_list <- function(x,...){
    if(class(x) == "otolith_contour"){
        r <- list(x)
        class(r) <- "otolith_contour_list"
        return(r)
    }
    if(class(x) != "list" & all(unlist(lapply(x,class)) == "otolith_image"))
        stop("This function can only be used on output from the read_image function or an otolith_contour object.")
    r <- do.call("c",x)
    class(r) <- "otolith_contour_list"
    return(r)
}

##' Extract otolith contours from images
##'
##' @param files file paths to read
##' @param ... Arguments passed to read_image
##' @return an otolith_contour_list object
##' @author Christoffer Moesgaard Albertsen
##' @seealso read_image
##' @export
read_images <- function(files, ...){
    ims <- lapply(files, read_image, ...)
    return(as.otolith_contour_list(ims))
}


##' @method c otolith_contour
##' @export
c.otolith_contour <- function(...){
    res <- list(...)
    class(res) <- "otolith_contour_list"
    return(res)
}

##' @method c otolith_contour_list
##' @export
c.otolith_contour_list <- function(...){
    res <- unlist(list(...), recursive = FALSE)
    class(res) <- "otolith_contour_list"
    return(res)
}

##' @method c otolith_image
##' @export
c.otolith_image <- function(...){
    lr1 <- lapply(list(...),as.otolith_contour_list)
    do.call("c",lr1)
}

