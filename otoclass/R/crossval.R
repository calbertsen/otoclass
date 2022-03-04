ConfusionMatrix <- function(G, P, accuracy = TRUE){
    ## Handle levels for factors!
    tab <- table(True = G, Predicted = P)
    if(!accuracy) return(tab)
    tab / rowSums(tab)[row(tab)]
}
metric_MatthewCorrelationCoefficient <- function(G, P){
    C <- ConfusionMatrix(G,P, FALSE)
    tk <- rowSums(C)
    pk <- colSums(C)
    cc <- sum(diag(C))
    s <- sum(C)
    (cc * s - sum(tk * pk)) / sqrt((s^2 - sum(pk^2)) * (s^2 - sum(tk^2)) )
}
metric_totalAccuracy <- function(G, P){
    mean(G==P)
}
metric_balancedAccuracy <- function(G, P){
    mean(diag(ConfusionMatrix(G,P)))
}


##' @export
leaveoneout <- function(f, ...)
    UseMethod("leaveoneout")

##' @export
leaveoneout.mlld <- function(f, prior = NULL){
    cl0 <- f$call
    cl0$group <- getGroupDataFrame(f)
    cl0$parlist <- f$pl
    cl0$doSdreport <- FALSE
    e <- environment(cl0)
    doOne <- function(i){
        cat(i,"\n")
        cl <- cl0
        cl$group[i,] <- NA
        suppressWarnings(oo <- capture.output(fi <- eval(cl,e)))
        predict(fi, prior = prior)$class[i]        
    }
    r <- sapply(seq_len(nrow(cl0$group)), doOne)
    tab <- table(True = cl0$group[,1], Predicted = r)    
    list(AccuracyMatrix = ConfusionMatrix(cl0$group[,1], r),
         MCC = metric_MatthewCorrelationCoefficient(cl0$group[,1], r),
         BalancedAccuracy = metric_balancedAccuracy(cl0$group[,1], r),
         TotalAccuracy = metric_totalAccuracy(cl0$group[,1], r),
         Predictions = r)
}



##' @export
crossval <- function(f,
                     folds,
                     ...){
    UseMethod("crossval")
}


## library(otoclass)
## group <- sample(1:2,1000, replace = TRUE)
## X <- t(sapply(group,function(i)c(-1,1)[i] + t(chol(matrix(c(1,0.8,0.8,1),2,2))) %*% rnorm(2)))
## x1 <- mlld(X,factor(group), featureSelection = "Lp", penalty = 2)
## x2 <- mlld(X,factor(group), featureSelection = "Lp", penalty = 2, estimateLambda=FALSE)
## cv1 <- crossval(x1)
## cv2 <- crossval(x2)

##' @export
crossval.mlld <- function(f,
                          folds,
                          prior = NULL,
                          ...){
    cl0 <- f$call
    e <- environment(cl0)
    cl0$group <- getGroupDataFrame(f)
    cl0$y <- eval(cl0$y,e)
    cl0$data <- eval(cl0$data,e)
    cl0$parlist <- f$pl
    cl0$doSdreport <- FALSE
    if(is.null(cl0$data)){
        cl0$data <- data.frame(ID = seq_len(length(cl0$group)))
    }
    if(!is.matrix(cl0$y))
        cl0$y <- matrix(cl0$y,ncol = 1)
    doOne <- function(i){
        cat(i,"\n")
        cl <- cl0
        cl$group <- cl0$group[-i,,drop=FALSE]
        cl$y <- cl0$y[-i,,drop=FALSE]
        cl$data <- cl0$data[-i,,drop=FALSE]
        suppressWarnings(oo <- capture.output(fi <- eval(cl,e)))
        predict(fi, y = cl0$y[i,,drop=FALSE], data = cl0$data[i,,drop=FALSE], prior = prior)$class        
    }
    kk <- sample(seq_len(folds), nrow(cl0$group), replace = TRUE)
    indx <- split(seq_len(nrow(cl0$group)), kk)
    r <- unlist(lapply(indx, doOne))[order(unlist(indx))]
    tab <- table(True = cl0$group[,1], Predicted = r)    
    list(AccuracyMatrix = ConfusionMatrix(cl0$group[,1], r),
         MCC = metric_MatthewCorrelationCoefficient(cl0$group[,1], r),
         BalancedAccuracy = metric_balancedAccuracy(cl0$group[,1], r),
         TotalAccuracy = metric_totalAccuracy(cl0$group[,1], r),
         Predictions = r,
         Fold = kk)
}


##' @export
forwardSelection <- function(f, ...)
    UseMethod("forwardSelection")


##' @export
forwardSelection.mlld <- function(f, prior = NULL, maxit = Inf, criterion = c("MCC","TA","BA"),testPct = 0.25){
    criterion <- match.arg(criterion)
    if(criterion == "MCC"){
        cFun <- metric_MatthewCorrelationCoefficient
    }else if(criterion == "TA"){
        cFun <- metric_totalAccuracy
    }else if(criterion == "BA"){
        cFun <- metric_balancedAccuracy
    }else{
        stop("Criteria not implemented")
    }
    cl0 <- f$call
    e <- environment(f$call)
    cl0$group <- getGroupDataFrame(f)
    cl0$y <- eval(cl0$y,e)
    cl0$data <- eval(cl0$data,e)
    if(!is.matrix(cl0$y) || (is.matrix(cl0$y) && ncol(cl0$y) <= 1))
        stop("Not enough features to select")
    cl0$doSdreport <- FALSE
    cl0$estimateLambda <- FALSE
    cl0$lp_penalty <- NA
    testIndex <- sort(unlist(lapply(split(seq_len(nrow(cl0$y)), cl0$group[,1]), function(ii) sample(ii,testPct * length(ii)))))
    trainIndex <- setdiff(seq_len(nrow(cl0$y)), testIndex)
    trY <- cl0$y[testIndex,,drop=FALSE]
    trD <- cl0$data[testIndex,,drop=FALSE]
    trG <- cl0$group[testIndex,,drop=FALSE]
    cl0$y <- cl0$y[trainIndex,,drop=FALSE]
    cl0$data <- cl0$data[trainIndex,,drop=FALSE]
    cl0$group <- cl0$group[trainIndex,,drop=FALSE]
    Selected <- c()
    Remain <- seq_len(ncol(cl0$y))
    Accuracy <- -Inf
    while(length(Remain) > 0 && length(Selected) <= maxit){
        cat("Selecting feature:",length(Selected)+1,"\n")
        ## Try all the models
        v <- sapply(Remain, function(inext){
            cl <- cl0
            cl$y <- cl0$y[, c(Selected, inext), drop=FALSE]
            suppressWarnings(oo <- capture.output(fi <- eval(cl,e)))
            P <- predict(fi,
                         y = trY[,c(Selected, inext),drop=FALSE],
                         data = trD,
                         prior = prior)$class
            cFun(trG[,1], P)
        })        
        v[is.nan(v) | is.na(v) | !is.finite(v)] <- -Inf
        ## Update
        if(any(v > tail(Accuracy,1))){
            cat("\t",sprintf("Criterion improved from %f to %f",tail(Accuracy,1),max(v)),"\n")
            cat("\t",sprintf("Selecting %s",f$muNames[[2]][Remain[which.max(v)]]),"\n")
            Selected <- c(Selected, Remain[which.max(v)])
            Remain <- setdiff(Remain,Selected)
            Accuracy <- c(Accuracy,max(v))
        }else{
            break;
        }
    }
    data.frame(Number = seq_along(Selected),
               Selected = f$muNames[[2]][Selected],
               Criterion = tail(Accuracy,-1),
               FeatureIndex = Selected)
}


##' @export
backwardSelection <- function(f, ...)
    UseMethod("backwardSelection")

##' @export
backwardSelection.mlld <- function(f, prior = NULL, maxit = Inf, criterion = c("MCC","TA","BA"),testPct = 0.25){
    criterion <- match.arg(criterion)
    if(criterion == "MCC"){
        cFun <- metric_MatthewCorrelationCoefficient
    }else if(criterion == "TA"){
        cFun <- metric_totalAccuracy
    }else if(criterion == "BA"){
        cFun <- metric_balancedAccuracy
    }else{
        stop("Criteria not implemented")
    }
    cl0 <- f$call
    e <- environment(f$call)
    cl0$group <- getGroupDataFrame(f)
    cl0$y <- eval(cl0$y,e)
    cl0$data <- eval(cl0$data,e)
    if(!is.matrix(cl0$y) || (is.matrix(cl0$y) && ncol(cl0$y) <= 1))
        stop("Not enough features to select")
    cl0$doSdreport <- FALSE
    cl0$estimateLambda <- FALSE
    cl0$lp_penalty <- NA
    testIndex <- sort(unlist(lapply(split(seq_len(nrow(cl0$y)), cl0$group[,1]), function(ii) sample(ii,testPct * length(ii)))))
    trainIndex <- setdiff(seq_len(nrow(cl0$y)), testIndex)
    trY <- cl0$y[testIndex,,drop=FALSE]
    trD <- cl0$data[testIndex,,drop=FALSE]
    trG <- cl0$group[testIndex,,drop=FALSE]
    cl0$y <- cl0$y[trainIndex,,drop=FALSE]
    cl0$data <- cl0$data[trainIndex,,drop=FALSE]
    cl0$group <- cl0$group[trainIndex,,drop=FALSE]
    Selected <- c()
    Remain <- seq_len(ncol(cl0$y))
    P0 <- predict(f,
                  y = trY,
                  data = trD,
                  prior = prior)$class    
    Accuracy <- cFun(trG[,1], P0)
    while(length(Remain) > 1 && length(Selected) <= maxit){
        cat("Removing feature:",length(Remain)+1,"\n")
        ## Try all the models
        v <- sapply(Remain, function(inext){
            cl <- cl0
            cl$y <- cl0$y[, setdiff(Remain, inext), drop=FALSE]
            suppressWarnings(oo <- capture.output(fi <- eval(cl,e)))
            P <- predict(fi,
                         y = trY[,setdiff(Remain, inext),drop=FALSE],
                         data = trD,
                         prior = prior)$class
            cFun(trG[,1], P)
        })
        v[is.nan(v) | is.na(v) | !is.finite(v)] <- -Inf
        ## Update
        if(any(v > tail(Accuracy,1))){
            cat("\t",sprintf("Criterion improved from %f to %f",tail(Accuracy,1),max(v)),"\n")
            cat("\t",sprintf("Removing %s",f$muNames[[2]][Remain[which.max(v)]]),"\n")
            Selected <- c(Selected, Remain[which.max(v)])
            Remain <- setdiff(Remain,Selected)
            Accuracy <- c(Accuracy,max(v))
        }else{
            break;
        }
    }
    list(Remaining = f$muNames[[2]][Remain],
         Removed = data.frame(Number = seq_along(Selected),
                              Selected = f$muNames[[2]][Selected],
                              Criterion = tail(Accuracy,-1),
                              FeatureIndex = Selected)
         )
}
