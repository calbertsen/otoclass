##' Calculate Elliptical Fourier Descriptors
##'
##' @param dat 
##' @param N 
##' @param n 
##' @param returnAsList 
##' @return List of elliptical fourier descriptors
##' @author Christoffer Moesgaard Albertsen
##' @export
efd <- function(dat,N,n=nrow(dat),returnAsList=FALSE){
    t <- seq(0,1,len=n)
    T <- t[length(t)]
    dx <- diff(dat[,1])
    dy <- diff(dat[,2])
    dt <- diff(t)
    a <- Vectorize(function(n){
        T/(2*n^2*pi^2)*sum(dx/dt*(cos(2*n*pi*t[-1]/T)-cos(2*n*pi*t[-length(t)]/T)))
    })
    b <-  Vectorize(function(n){
        T/(2*n^2*pi^2)*sum(dx/dt*(sin(2*n*pi*t[-1]/T)-sin(2*n*pi*t[-length(t)]/T)))
    })
    cc <- Vectorize(function(n){
        T/(2*n^2*pi^2)*sum(dy/dt*(cos(2*n*pi*t[-1]/T)-cos(2*n*pi*t[-length(t)]/T)))
    })
    d <- Vectorize(function(n){
        T/(2*n^2*pi^2)*sum(dy/dt*(sin(2*n*pi*t[-1]/T)-sin(2*n*pi*t[-length(t)]/T)))
    })

    acalc <- a(1:N)
    bcalc <- b(1:N)
    ccalc <- cc(1:N)
    dcalc <- d(1:N)
    a0 <- 1/T * sum(dat[,1]*c(0,dt))
    names(a0) <- "A0"
    c0 <- 1/T * sum(dat[,2]*c(0,dt))
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

##' Calculate XY-coordinates from Elliptical Fourier Descriptors
##'
##' @param n 
##' @param A 
##' @param B 
##' @param C 
##' @param D 
##' @param A0 
##' @param C0 
##' @return Matrix of coordinates
##' @author Christoffer Moesgaard Albertsen
##' @export
efd2coord <- function(n,A,B,C,D,A0=0,C0=0){
    tt <- seq(0,1,len=n)
    T <- tt[length(tt)]
    x <- A0 + sapply(tt,function(t)sum(A*cos(2*1:length(A)*pi*t/T)+B*sin(2*1:length(A)*pi*t/T)))
    y <- C0 + sapply(tt,function(t)sum(C*cos(2*1:length(A)*pi*t/T)+D*sin(2*1:length(A)*pi*t/T)))
    return(cbind(x,y))
}
