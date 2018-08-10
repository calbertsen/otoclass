

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
##' @param nCuts Function that calculates number of bins used for calculating correlations
##' @return A vector of selected features
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stats na.omit
##' @export
fcbf <- function(train,
                  group,
                  delta,
                  lambda=1,
                  nCuts = function(x){round(length(unique(x)) / 5)}){
    # http://www.cs.waikato.ac.nz/ml/publications/1997/Hall-LSmith97.pdf
                                        # http://pdf.aminer.org/000/335/746/feature_selection_for_high_dimensional_data_a_fast_correlation_based.pdf
    ## Exclude variables with no variance
    Slist <- numeric(0)
    vals <- apply(train,2,function(x){
        SU(x,group,nx = max(2,nCuts(x)))
    })
    Slist <- order(vals,decreasing=TRUE)
    ## Exclude variables with no variance
    hasVar <- apply(train,2,function(x)!isTRUE(all.equal(sd(x),0)))
    Slist <- Slist[vals[Slist]>=delta & hasVar[Slist]]
    if(length(Slist) == 0)
        stop("There are no variables that are more than delta correlated to the groups.")
    if(length(Slist) == 1)
        return(Slist)
    N <- length(Slist)
    for(p in 1:N){
        if(!is.na(Slist[p])){
            indxUse <- na.omit(Slist[-(1:p)])
            trainUse <- train[,indxUse]
            if(class(trainUse)=="numeric"){
                nT1 <- max(2,nCuts(trainUse))
                nT2 <- max(2,nCuts(train[,Slist[p]]))
                valp <- SU(trainUse,train[,Slist[p]],nT1,nT2)
                Slist[indxUse[valp >= lambda*vals[stats::na.omit(Slist[-(1:p)])]]] <- NA
            }else if(dim(trainUse)[2]>0){
                valp <- apply(train[,stats::na.omit(Slist[-(1:p)])],2,function(x){
                    nT1 <- max(2,nCuts(x))
                    nT2 <- max(2,nCuts(train[,Slist[p]]))
                    SU(x,train[,Slist[p]],nT1,nT2)
                })
                Slist[Slist%in%indxUse[valp >= lambda*vals[stats::na.omit(Slist[-(1:p)])]]] <- NA
            }
        }        
    }
     return(as.vector(na.omit(Slist)))
}

