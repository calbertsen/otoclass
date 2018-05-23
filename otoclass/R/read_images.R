
##' @export
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

##' @export
read_images <- function(files, ...){
    ims <- lapply(files, read_image, ...)
    return(as.otolith_contour_list(ims))
}
