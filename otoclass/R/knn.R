
#' k Nearest Neighbour
#'
#' \code{knn} returns the votes to each level factor(group) from the kn nearest neighbours to test from train.
#' @param train A matrix of (continuous) features for the training sample
#' @param group A vector of group labels for the training sample
#' @param test A matrix of (continuous) features for the test sample
#' @param kn Number of neighbours to consider
#' @param dist Distance measure to use ("L1", "L2","Linf"). Default is Euclidian distance ("L2").

#' @export
knn <- function(train,
                group,
                test,
                kn = 1,
                dist = c("L2", "L1","Linf")
                ) {
    if(!is.factor(group)){
        group <- factor(group)
    }
    dist <- match.arg(dist)
    res <- list()
    class(res) <- "oto_knn"
    res$probabilities <- .Call('knn', PACKAGE = 'otoclass', train = t(train), group = group, test = t(test), kn = as.integer(kn), disttype = as.integer(factor(dist,c("L2", "L1","Linf"))))
    colnames(res$probabilities) <- levels(group)
    rownames(res$probabilities) <- rownames(test)
    res$predicted <- factor(levels(group)[apply(res$probabilities,1,which.max)], levels = levels(group))

    return(res)
}
