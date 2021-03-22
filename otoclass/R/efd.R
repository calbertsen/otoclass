##' Calculate Elliptical Fourier Descriptors
##'
##' @param x contour object 
##' @param N Number of harmonics to calculate
##' @param returnAsList Should result be returned as a list that can be passed to \code{efd2coord}?
##' @param normalize Normalize using method from Ferson et al. 1985?
##' @param ... other argiments
##' @return List of elliptical fourier descriptors
##' @author Christoffer Moesgaard Albertsen
##' @references Kuhl, F. P., Giardina, C. R. (1981) Elliptic Fourier Features of a Closed Contour. Computer Graphics and Image Processing. 18:236-258
##' Ferson, S., Rohlf, F. J., Koehn, R. K. (1985) Measuring Shape Variation of Two-Dimensional Outlines. Systematic Zoology. 34(1):59-68
##' @examples
##' Ellipsis <- function(theta,r1=1,r2=2){
##'   xx <- seq(0,2*pi,len=100)
##'   X <- cbind(r1*cos(xx),r2*sin(xx))
##'   R <- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
##'   Y <- t(apply(X,1,function(x) as.vector(R %*% x)))
##'   return(Y)
##' }
##' efd(Ellipsis(pi/2),1)
##'
##' efd(Ellipsis(pi/2),1,normalize = TRUE)
##'
##' efd(Ellipsis(-pi/2),1,normalize = TRUE)
##'
##' @export
efd <- function(x, N, ...){
    UseMethod("efd")
}

##' @rdname efd
##' @method efd matrix
##' @export
efd.matrix <- function(x, N, returnAsList = FALSE, normalize = FALSE, ...){
    if(!isTRUE(all.equal(x[1,],x[nrow(x),])))
        x <- rbind(x,x[1,])
    xu <- x[c(TRUE,diff(x[,1]) != 0 ) | c(TRUE,diff(x[,2]) != 0 ), ]
    res <- .Call("efd",
                 dat = xu,
                 N = as.integer(N),
                 normalize = normalize,
                 PACKAGE="otoclass")
    names(res$A0) <- "A0"
    names(res$C0) <- "C0"
    names(res$A) <- paste0("A",1:N)
    names(res$B) <- paste0("B",1:N)
    names(res$C) <- paste0("C",1:N)
    names(res$D) <- paste0("D",1:N)

    if(returnAsList){
        return(res)
    }else{
        names(res) <- NULL
        return(unlist(res))
    }
}


##' @rdname efd
##' @method efd otolith_contour
##' @export
efd.otolith_contour <- function(x, N, ...){
    efd.matrix(x = unclass(x), N = N, ...)
}

##' @rdname efd
##' @method efd otolith_image
##' @export
efd.otolith_image <- function(x, N, ...){
    lapply(x, efd, N=N, ...)
}

##' @rdname efd
##' @method efd otolith_contour_list
##' @export
efd.otolith_contour_list <- function(x, N, ...){
    lapply(x, efd, N=N, ...)    
}

##' Calculate XY-coordinates from Elliptical Fourier Descriptors
##'
##' @param n Number of coordinates to calculate
##' @param A A EFD coefficients
##' @param B B EFD coefficients
##' @param C C EFD coefficients
##' @param D D EFD coefficients
##' @param A0 A0 EFD coefficients
##' @param C0 C0 EFD coefficients
##' @return Matrix of coordinates
##' @author Christoffer Moesgaard Albertsen
efd2coordOld <- function(n,A,B,C,D,A0=0,C0=0){
    tt <- seq(0,1,len=n)
    T <- tt[length(tt)]
    x <- A0 + sapply(tt,function(t)sum(A*cos(2*1:length(A)*pi*t/T)+B*sin(2*1:length(A)*pi*t/T)))
    y <- C0 + sapply(tt,function(t)sum(C*cos(2*1:length(A)*pi*t/T)+D*sin(2*1:length(A)*pi*t/T)))
    return(cbind(x,y))
}


##' Calculate XY-coordinates from Elliptical Fourier Descriptors
##'
##' @param n Number of coordinates to calculate
##' @param A A EFD coefficients
##' @param B B EFD coefficients
##' @param C C EFD coefficients
##' @param D D EFD coefficients
##' @param A0 A0 EFD coefficients
##' @param C0 C0 EFD coefficients
##' @return Matrix of coordinates
##' @author Christoffer Moesgaard Albertsen
##' @export
efd2coord <- function(n,A,B,C,D,A0=0,C0=0){
    coords <- .Call("efd2coordSEXP",
                    efd = cbind(A,B,C,D),
                    N=as.integer(n),
                    A0 = A0,
                    C0 = C0,
                    PACKAGE = "otoclass")
    return(coords)
}


 Ellipsis <- function(theta,r1=1,r2=2,len=100, x0=0,y0=0){
   xx <- seq(0,2*pi,len=len)
   X <- cbind(r1*cos(xx),r2*sin(xx))
   R <- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
   Y <- t(apply(X,1,function(x) as.vector(R %*% x))) + matrix(c(x0,y0),len,2, byrow=TRUE)
   return(Y)
 }

## xx <- Ellipsis(pi/2) + matrix(c(-2,5),100,2,byrow=TRUE)
## .Call("efd",xx,1L,FALSE,PACKAGE="otoclass")

efdOld <- function(dat, N, returnAsList=FALSE, normalize = FALSE){
    if(!all(dat[1,] == dat[nrow(dat),]))
        dat <- rbind(dat,dat[1,])
    n <- nrow(dat)
 
    dx <- diff(dat[,1])
    dy <- diff(dat[,2])
    dt <- sqrt(dx^2 + dy^2)
    tt <- c(0,cumsum(dt))
    T <- tail(tt,1)
    a <- Vectorize(function(n){
        T/(2*n^2*pi^2)*sum(dx/dt*(cos(2*n*pi*tt[-1]/T)-cos(2*n*pi*tt[-length(tt)]/T)))
    })
    b <-  Vectorize(function(n){
        T/(2*n^2*pi^2)*sum(dx/dt*(sin(2*n*pi*tt[-1]/T)-sin(2*n*pi*tt[-length(tt)]/T)))
    })
    cc <- Vectorize(function(n){
        T/(2*n^2*pi^2)*sum(dy/dt*(cos(2*n*pi*tt[-1]/T)-cos(2*n*pi*tt[-length(tt)]/T)))
    })
    d <- Vectorize(function(n){
        T/(2*n^2*pi^2)*sum(dy/dt*(sin(2*n*pi*tt[-1]/T)-sin(2*n*pi*tt[-length(tt)]/T)))
    })

    acalc <- a(1:N)
    bcalc <- b(1:N)
    ccalc <- cc(1:N)
    dcalc <- d(1:N)
    a0 <- 1/T * sum(dat[,1]*c(0,dt))
    c0 <- 1/T * sum(dat[,2]*c(0,dt))

    if(normalize){
        theta <- (0.5 * atan2( 2 * ( acalc[1] * bcalc[1] + ccalc[1] * dcalc[1] ),
        (acalc[1]^2 - bcalc[1]^2 + ccalc[1]^2 - dcalc[1]^2)))
        as <- acalc[1] * cos(theta) + bcalc[1] * sin(theta)
        cs <- ccalc[1] * cos(theta) + dcalc[1] * sin(theta)
        scale <-  1 / sqrt(as^2 + cs^2)
        phi <- atan2(cs, as)
        T1 <- matrix(c(cos(phi),-sin(phi),sin(phi),cos(phi)),2,2)
        T2 <- function(n) matrix(c(cos(n * theta),sin(n * theta),-sin(n * theta),cos(n * theta)),2,2)
        nefd <- do.call("rbind",sapply(1:N, function(i){ as.vector(t(scale * T1 %*% matrix(c(acalc[i],ccalc[i],bcalc[i],dcalc[i]),2,2) %*% T2(i)))},simplify = FALSE))
        acalc <- nefd[,1]
        bcalc <- nefd[,2]
        ccalc <- nefd[,3]
        dcalc <- nefd[,4]
        a0 <- 0
        c0 <- 0
    }
    
    names(a0) <- "A0"
    names(c0) <- "C0"
    names(acalc) <- paste0("A",1:N)
    names(bcalc) <- paste0("B",1:N)
    names(ccalc) <- paste0("C",1:N)
    names(dcalc) <- paste0("D",1:N)
    res <- list(A = acalc,
                B = bcalc,
                C = ccalc,
                D = dcalc,
                A0 = a0,
                C0 = c0)

    if(returnAsList){
        return(res)
    }else{
        names(res) <- NULL
        return(unlist(res))
    }
}


##' Fourier Series Basis for Splines
##'
##' @param x the predictor variable
##' @param N Half the number of basis functions. Both a cos and a sin basis function will be added.
##' @param Boundary Boundary points that determines a period
##' @return A matrix of dimension 'c(length(x), 2*N)'.
##' @author Christoffer Moesgaard Albertsen
##' @export
fs <- function(x, N, Boundary = range(x)){
    if(length(N) == 1)
        N <- seq_len(N)
    Boundary <- sort(Boundary)
    x0 <- (x - Boundary[1]) / diff(Boundary)
    basis <- do.call("cbind",lapply(N, function(i){
        m <- cbind(A = cos(2 * pi * i * x0),
                   B = sin(2 * pi * i * x0))
        colnames(m) <- paste0(colnames(m),i)
        m
    }))
    a <- list(N = N, Boundary = Boundary)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("fs","matrix")
    basis
}
##' @author Christoffer Moesgaard Albertsen
##' @method makepredictcall fs
##' @importFrom stats makepredictcall
##' @export
makepredictcall.fs <- function(var, call){
    if (as.character(call)[1L] != "fs") 
        return(call)
     at <- attributes(var)[c("N", "Boundary")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx

}
