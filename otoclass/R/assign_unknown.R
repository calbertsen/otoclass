##' Assign unknown contours by binomial regression on elliptic fourier descriptors
##'
##' @param nimList A list of normalized contours
##' @param N Number of harmonics for the fourier descriptors
##' @param delta passed to fcbf for feature selection
##' @param lambda passed to fcbf for feature selection
##' @return A list of normalized contours where 'Unknown' Position attributes are replaced by 'Left' or 'Right'.
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats glm predict.glm binomial
##' @export
assign_unknown <- function(nimList,N=60,delta=0.01,lambda=1, link = "probit"){
    efds <- do.call("rbind",lapply(nimList,efd,N=60,returnAsList=FALSE))
    position <- unlist(lapply(nimList,attr,which="Position"))
    indx <- which(position == "Unknown")
    featUse <- fcbf(efds[-indx,],factor(position[-indx]),delta,lambda)
    pgrp <- factor(position[-indx])
    fit <- stats::glm(pgrp ~ .,
               family=stats::binomial(link),
               data=data.frame(efds[-indx,featUse]))
    pred <- predict(fit,newdata=data.frame(efds[indx,featUse]),type="response")
    for(i in 1:length(indx)){
        attr(nimList[[indx[i]]],"Position") <- levels(pgrp)[as.numeric(pred[i] > 0.5) + 1]
        attr(nimList[[indx[i]]],"Position_probability_Right") <- pred[i]
    }
    attr(nimList,"Assign_unknown_Features_used") <- colnames(efds)[featUse]
    attr(nimList,"Assign_unknown_N") <- N
    attr(nimList,"Assign_unknown_Delta") <- delta
    attr(nimList,"Assign_unknown_Lambda") <- lambda
    attr(nimList,"Assign_unknown_Link") <- link
    return(nimList)
}