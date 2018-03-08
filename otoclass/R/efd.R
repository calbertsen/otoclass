##' Calculate Elliptical Fourier Descriptors
##'
##' @param dat 
##' @param N 
##' @param n 
##' @param returnAsList 
##' @param normalize Normalize using method from Ferson et al. 1985?
##' @return List of elliptical fourier descriptors
##' @author Christoffer Moesgaard Albertsen
##' @export
efd <- function(dat,N,n=nrow(dat),returnAsList=FALSE, normalize = FALSE){
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

    if(normalize){
        theta <- (0.5 * atan((2 * (acalc[1] * bcalc[1] + ccalc[1] * dcalc[1] )) /
                          (acalc[1]^2 + bcalc[1]^2 + ccalc[1]^2 + dcalc[1]^2))) %% pi
        as <- acalc[1] * cos(theta) + bcalc[1] * sin(theta)
        cs <- ccalc[1] * cos(theta) + dcalc[1] * sin(theta)
        scale <-  1 / sqrt(as^2 + cs^2)
        phi <- atan(as / cs) %% pi
        T1 <- matrix(c(cos(phi),-sin(phi),sin(phi),cos(phi)),2,2)
        T2 <- function(n) matrix(c(cos(n * theta),-sin(n * theta),sin(n * theta),cos(n * theta)),2,2)
        nefd <- do.call("rbind",sapply(1:N, function(i){ as.vector(t(scale * T1 %*% matrix(c(acalc[i],ccalc[i],bcalc[i],dcalc[i]),2,2) %*% T2(i)))},simplify = FALSE))
        acalc <- nefd[,1]
        bcalc <- nefd[,2]
        ccalc <- nefd[,3]
        dcalc <- nefd[,4]
    }
    
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
