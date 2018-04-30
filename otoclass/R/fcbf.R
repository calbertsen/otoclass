

SU <- function(x,y,nx=5,ny=5){
    if(is.factor(x)){
        xgrp <- x
    }else{
        xgrp <- cut(x,nx)
    }
    if(is.factor(y)){
        ygrp <- y
    }else{
        ygrp <- cut(y,ny)
    }
    px <- table(xgrp)/length(xgrp)
    py <- table(ygrp)/length(ygrp)
    pxy <- apply(table(xgrp,ygrp),2,function(x)x/sum(x)) 
    lpx <- ifelse(px==0,0,log(px,base=2))
    lpy <- ifelse(py==0,0,log(py,base=2))
    lpxy <- ifelse(pxy==0,0,log(pxy,base=2))
    Hx <- -sum(px*lpx)
    Hy <- -sum(py*lpy)
    Hxy <- -sum(py*colSums(pxy*lpxy))
    return(2*(Hx-Hxy)/(Hx+Hy))
}

##' Fast Correlation Based Filter for Feature Selection in High-Dimensional Data
##'
##' @references http://pdf.aminer.org/000/335/746/feature_selection_for_high_dimensional_data_a_fast_correlation_based.pdf
##' @param train A matrix of training data
##' @param group A vector of groups
##' @param delta Tuning parameter
##' @param lambda Tuning parameter
##' @return A vector of selected features
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats na.omit
##' @export
fcbf <- function(train,group,delta,lambda=1){
    # http://www.cs.waikato.ac.nz/ml/publications/1997/Hall-LSmith97.pdf
    # http://pdf.aminer.org/000/335/746/feature_selection_for_high_dimensional_data_a_fast_correlation_based.pdf
    Slist <- numeric(0)
    vals <- apply(train,2,function(x)SU(x,group))
    Slist <- order(vals,decreasing=TRUE)
    Slist <- Slist[vals[Slist]>=delta]
    N <- length(Slist)
    for(p in 1:N){
        if(!is.na(Slist[p])){
            indxUse <- na.omit(Slist[-(1:p)])
            trainUse <- train[,indxUse]
            if(class(trainUse)=="numeric"){
                valp <- SU(trainUse,train[,Slist[p]])
                Slist[indxUse[valp >= lambda*vals[stats::na.omit(Slist[-(1:p)])]]] <- NA
            }else if(dim(trainUse)[2]>0){
                valp <- apply(train[,stats::na.omit(Slist[-(1:p)])],2,function(x)SU(x,train[,Slist[p]]))
                Slist[Slist%in%indxUse[valp >= lambda*vals[stats::na.omit(Slist[-(1:p)])]]] <- NA
            }
        }        
    }
    return(as.vector(na.omit(Slist)))
}

