
calculateVarBetweenGroups <- function(x, ...){
    UseMethod("calculateVarBetweenGroups")
}

##' @export
calculateVarBetweenGroups.mlld <- function(x, ...){
    if(!is.na(x$varBetweenGroups))
        return(x$varBetweenGroups)
    X <- model.matrix(x$terms,data = x$data)
    beta <- x$pl$mu
    group <- x$group
    muList <- sapply(1:nrow(X), function(i){
        X[i,] %*% beta[,,group[i]]
    }, simplify = FALSE)
    muAll <- colMeans(x$train)
    Sigma <- Reduce("+",lapply(muList, function(xx) t(xx-muAll) %*% (xx-muAll))) / (nrow(X)-1)
    return(Sigma)    
}


calculateVarWithinGroups <- function(x, ...){
    UseMethod("calculateVarWithinGroups")
}

##' @export
calculateVarWithinGroups.mlld <- function(x, ...){
    if(!is.na(x$varWithinGroups))
        return(x$varWithinGroups)
    sigmas <- exp(x$pl$logSigma)
    corpars <- x$pl$corpar
    Sigma <- Reduce("+",sapply(1:ncol(sigmas), function(i){
        A <- diag(sigmas[,i], nrow(sigmas), nrow(sigmas))
        M <- diag(1,nrow(A),nrow(A))
        M[lower.tri(M)] <- corpars[,i]
        A%*%(M%*%t(M))%*%t(A)
    }, simplify = FALSE)) / ncol(sigmas)
    return(Sigma)
}




##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param x 
##' @param ... 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @export
projection <- function(x, ...){
    UseMethod("projection")
}

##' @export
projection.mlld <- function(x, test, ...){
    if(missing(test)){
        test <- x$train
        group <- x$group
    }else{
        group <- factor(rep(NA,nrow(test)),levels = levels(x$group))
    }
    ee <- eigen(solve(x$varWithinGroups) %*% x$varBetweenGroups,symmetric = FALSE)
    noImPart <- sapply(ee$values, function(x) isTRUE(all.equal(Im(x),0)))
    nonZeroRe <- !sapply(ee$values / sum(ee$values),function(x)isTRUE(all.equal(Re(x),0)))
    eeUse <- nonZeroRe & noImPart
    if(any(!noImPart & nonZeroRe))
        warning("Imaginary part of eigen vectors ignored.")
    evec <- Re(ee$vectors[,eeUse])
    if(is.data.frame(test))
        test <- as.matrix(test)
    projdat <- test %*% evec
    colnames(projdat) <- paste0("CD",1:ncol(projdat))
    rownames(projdat) <- rownames(test)
    attr(projdat,"relative_importance") <- Re(ee$values[eeUse]) / sum(Re(ee$values[eeUse]))
    attr(projdat,"group") <- group
    class(projdat) <- "projection"
    return(projdat)
}



addTrans <- Vectorize(function(name,alpha=1){
    if(is.na(name) | is.null(name))
        return(name)
    arg <- as.list(grDevices::col2rgb(name)/255)
    names(arg) <- c("red","green","blue")
    arg$alpha <- alpha
    do.call(grDevices::rgb,arg)
})


##' @export
plot.projection <- function(x, onlyFirst = FALSE, ...){
    group <- attr(x,"group")
    lvl <- levels(attr(x,"group"))
    if(any(is.na(group))){
        group <- as.character(group)
        group[is.na(group)] <- "Unknown"
        group <- factor(group, levels = c(lvl,"Unknown"))
    }
    args <- list(...)
    if(!is.null(args$col)){
        args$col <- rep(args$col,nlevels(group))
    }else{
        args$col <- addTrans(1:nlevels(group),0.3)
    }
    if(ncol(x) == 1 | onlyFirst){
        pdList <- split(x, group)
        dens <- lapply(pdList,density, cut = 0)
        xrng <- range(unlist(lapply(dens,function(xx)xx$x)))
        yrng <- range(unlist(lapply(dens,function(xx)xx$y)))
        if(is.null(args$xlim))
            args$xlim <- xrng
        if(is.null(args$ylim))
            args$ylim <- yrng
        plot(0,0, type = "n",
             xlab = sprintf("CD1 (%.1f %%)",attr(x,"relative_importance")[1]*100),
             ylab = "Density",
             ...)
        for(i in 1:length(dens)){
            arg0 <- args
            arg0$x <- dens[[i]]
            arg0$col <- args$col[i]
            do.call(graphics::lines, arg0)
        }
    }else{
        args$col <- args$col[group]
        if(is.null(args$pch))
            args$pch <- 16
        if(is.null(args$xlab))
            args$xlab <- sprintf("CD1 (%.1f %%)",attr(x,"relative_importance")[1]*100)
        if(is.null(args$ylab))
            args$ylab <- sprintf("CD2 (%.1f %%)",attr(x,"relative_importance")[2]*100)
        args$x <- unclass(x)
        do.call(graphics::plot,args)
    }
}
