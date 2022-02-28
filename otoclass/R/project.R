


calculateVarBetweenGroups <- function(x, ...){
    UseMethod("calculateVarBetweenGroups")
}

##' @export
calculateVarBetweenGroups.mlld <- function(x, ...){
    if(!is.na(x$varBetweenGroups))
        return(x$varBetweenGroups)
    grp <- getGroupDataFrame.mlld(x)
    pr <- predict(x, group = grp)
    muList <- lapply(lapply(seq_along(levels(pr$class)), function(i) {
        vv <- pr$fit[as.numeric(pr$class)==i,,i, drop = TRUE]
        matrix(vv,sum(as.numeric(pr$class)==i),dim(pr$fit)[2])
    }), colMeans)
    muAll <- colMeans(x$y)
    N <- length(muList)
    Sigma <- Reduce("+",lapply(muList, function(xx) (xx-muAll) %*% t(xx-muAll))) / N
    colnames(Sigma) <- rownames(Sigma) <- names(muAll)
    return(Sigma)    
}


calculateVarWithinGroups <- function(x, ...){
    UseMethod("calculateVarWithinGroups")
}

##' @export
calculateVarWithinGroups.mlld <- function(x, ...){
    if(!is.na(x$varWithinGroups))
        return(x$varWithinGroups)
    vcx <- vcov(x)
    Sigma <- Reduce("+",vcx) / length(vcx)
    return(Sigma)
}




##' LDA-like projection
##'
##' @param x 
##' @param ... 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @export
projection <- function(x, ...){
    UseMethod("projection")
}

##' @export
projection.mlld <- function(x, y, ...){
    if(missing(y))
        y <- x$y
        ## group <- x$group
    ## }else if(!(missing(y) & missing(group))){
    ##     stop("Both y and group must be given.")
    ## }
    ## if(ncol(y) != ncol(x$y)){
    ##     stop("The new y does not match the fitted object.")
    ## }
    ## if(!isTRUE(all.equal(levels(group), levels(x$group)))){
    ##     stop("The new group does not match the fitted object.")
    ## }
    ## if(nrow(y) != length(group))
    ##     stop("y and group does not match.")
    ee <- eigen(solve(x$varWithinGroups) %*% x$varBetweenGroups,symmetric = FALSE)
    noImPart <- sapply(ee$values, function(x) isTRUE(all.equal(Im(x),0)))
    nonZeroRe <- !sapply(ee$values / sum(ee$values),function(x)isTRUE(all.equal(Re(x),0)))
    if(any(!noImPart & nonZeroRe))
        warning("Imaginary part of eigen vectors ignored.")
    eeUse <- head(which(nonZeroRe & noImPart), pmin(length(x$muNames[[3]]),ncol(y))-1)
    evec <- Re(ee$vectors[,eeUse, drop = FALSE])
    rownames(evec) <- x$muNames[[2]]
    colnames(evec) <- paste0("CD",1:ncol(evec))
    projdat <- y %*% evec
    colnames(projdat) <- paste0("CD",1:ncol(projdat))
    rownames(projdat) <- rownames(y)
    attr(projdat,"relative_importance") <- Re(ee$values[eeUse]) / sum(Re(ee$values[eeUse]))
    ## attr(projdat,"group") <- group
    attr(projdat,"directions") <- evec
    attr(projdat,"group") <- getGroupDataFrame(x)
    class(projdat) <- "projection"    
    return(projdat)
}

##' @export
print.projection <- function(x){
    cat("Direction\n")
    print(attr(x,"directions"))
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
plot.projection <- function(x, onlyFirst = FALSE, xlim = xrng, ylim = yrng, ...){
    group <- attr(x,"group")[,1]
    lvl <- levels(attr(x,"group"))[,1]
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
    if(!is.null(args$border)){
        args$border <- rep(args$border,nlevels(group))
    }else{
        args$border <- rep(NA,nlevels(group))
    }
    if(ncol(x) == 1 | onlyFirst){
        pdList <- split(x, group)
        dens <- lapply(pdList,density)
        xrng <- range(unlist(lapply(dens,function(xx)xx$x)))
        yrng <- range(unlist(lapply(dens,function(xx)xx$y)))
        if(is.null(xlim))
            xrng <- xlim
        if(is.null(ylim))
            yrng <- ylim
        plot(0,0, type = "n",
             xlab = sprintf("CD1 (%.1f %%)",attr(x,"relative_importance")[1]*100),
             ylab = "Density",
             ylim = yrng,
             xlim = xrng,             
             ...)
        for(i in 1:length(dens)){
            arg0 <- args
            arg0$x <- c(dens[[i]]$x,rev(dens[[i]]$x))
            arg0$y <- c(dens[[i]]$y,rev(dens[[i]]$y * 0))
            arg0$col <- args$col[i]
            arg0$border <- args$border[i]
            do.call(graphics::polygon, arg0)
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
