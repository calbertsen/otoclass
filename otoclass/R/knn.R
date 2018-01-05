
#' k Nearest Neighbour
#'
#' \code{knn} returns the votes to each level factor(group) from the kn nearest neighbours to test from train.
#' @param train A matrix of (continuous) features for the training sample
#' @param group A vector of group labels for the training sample
#' @param test A matrix of (continuous) features for the test sample
#' @param kn Number of neighbours to consider
#' @param dist Distance measure to use ("euclidian", "manhattan","chebyshev")

#' @export
knn <- function(train, group, test, kn, dist = "euclidian") {
    if(!is.factor(group)){
        group <- factor(group)
    }
    disttype <- as.numeric(factor(dist,levels=c("euclidian", "manhattan","L infinity")))
    if(is.na(disttype)) stop("Invalid distance measure. Must be one of: euclidian, manhattan, and chebyshev")
    res <- list()
    class(res) <- "oto_knn"
    res$probabilities <- .Call('otoclass_knn', PACKAGE = 'otoclass', train, as.numeric(group)-1, test, kn, disttype)
    colnames(res$probabilities) <- levels(group)
    rownames(res$probabilities) <- rownames(test)
    res$predicted <- factor(levels(group)[apply(res$probabilities,1,which.max)], levels = levels(group))

    return(res)
}
