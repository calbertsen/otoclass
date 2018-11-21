##' @export
crossval <- function(x,
                     folds,
                     by,
                     type = c("confusion", "feature"),
                     nPerBy = 5,
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
crossval.mlld <- function(x,
                          folds = nrow(eval(x$call$train)),
                          by,
                          type = c("confusion", "feature"),
                          nPerBy = 5,
                          gridN = 20,
                          gridMaxStep = 5,
                          ...){
    type <- match.arg(type)
    if(type != "confusion")
        stop("Cross validation type not implemented yet!")
    ## Sample 
    if(missing(by)){
        N <- nrow(eval(x$call$train))
        foldIndex <- sample(1:folds,size = N, replace = N > folds)
    }else{
        stop("Stratified sampling not implemented yet!")
    }

    if(type == "confusion"){
        cl <- as.list(x$call)
        pl <- x$pl
        trueGroup <- eval(cl$group)
        res <- sapply(sort(unique(foldIndex)), function(i){
            FUN <- eval(cl[[1]])
            args <- lapply(cl[-1],function(aa)eval(aa))
            args$test <- args$train[foldIndex == i, , drop = FALSE]
            args$train <- args$train[foldIndex != i, , drop = FALSE]
            args$group <- args$group[foldIndex != i]
            args$parlist <- pl
            if(!is.null(args$data)){
                args$dataTest <- args$data[foldIndex == i, , drop = FALSE]
                args$data <- args$data[foldIndex != i, , drop = FALSE]
            }
            args$estimateUnbiasedTestProportions <- FALSE
            args$biasCorrectionGroup <- NULL
            suppressWarnings(y <- do.call(FUN,args))
            pred <- predict(y)$class
            cbind(data.frame(true=trueGroup[foldIndex == i], prediction=pred),args$dataTest)
        }, simplify = FALSE)
        res <- do.call("rbind",res)
        op <- list(prediction = res,
                   crossval_table = table(True=res$true,Predicted=res$prediction))
        op$total_success <- sum(diag(op$crossval_table)) / sum(op$crossval_table)
        op$groupwise_success <- diag(op$crossval_table) / rowSums(op$crossval_table)
        op$group_numbers <- rowSums(op$crossval_table)
        op$confusion_matrix <- op$crossval_table / rowSums(op$crossval_table)
        class(op) <- c("crossval_confusion","crossval_mlld")
        return(op)
    }else if(type == "feature"){
        warning("crossval for features should not be used yet")
        cl <- as.list(x$call)
        pl <- x$pl
        trueGroup <- eval(cl$group)
        res <- sapply(sort(unique(foldIndex)), function(i){
            FUN <- eval(cl[[1]])
            args <- lapply(cl[-1],function(aa)eval(aa))
            args$test <- args$train[foldIndex == i, , drop = FALSE]
            args$train <- args$train[foldIndex != i, , drop = FALSE]
            args$group <- args$group[foldIndex != i]
            args$parlist <- pl
            if(!is.null(args$data)){
                args$dataTest <- args$data[foldIndex == i, , drop = FALSE]
                args$data <- args$data[foldIndex != i, , drop = FALSE]
            }
            args$estimateLambda <- FALSE
            xvals <- seq(-gridMaxStep,gridMaxStep,len = gridN)
            sapply(xvals,function(ii){
                arg0 <- args
                arg0$parlist$logLambda <- arg0$parlist$logLambda + ii
                suppressWarnings(y <- do.call(FUN,args))
                pred <- predict(y)$class
                tab <- table(true=trueGroup[foldIndex == i], pred)
                min(diag(tab)/rowSums(tab))
            })
        }, simplify = FALSE)
        return(res)
    }
}
