
##' Print otolith image contour
##'
##' @param x Otolith image contour
##' @param ... further arguments
##' @return nothing
##' @author Christoffer Moesgaard Albertsen
##' @method print otolith_image
print.otolith_image <- function(x, ...){
    cat("Otolith contour of", attr(x,"File"),"\n")
    cat("Position:",attr(x,"Position"))
    invisible(x)
}
