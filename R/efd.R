#'
#' @export

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
    names(acalc) <- paste0("A",1:N)
    names(bcalc) <- paste0("B",1:N)
    names(ccalc) <- paste0("C",1:N)
    names(dcalc) <- paste0("D",1:N)
    res <- list(A = acalc,
                B = bcalc,
                C = ccalc,
                D = dcalc)

    if(returnAsList){
        return(res)
    }else{
        names(res) <- NULL
        return(unlist(res))
    }
}

#' @export

efd2coord <- function(n,a,b,cc,d){
    t <- seq(0,1,len=n)
    T <- t[length(t)]
    x <- sum(a*cos(2*1:length(a)*pi*t/T)+b*sin(2*1:length(a)*pi*t/T))
    y <- sum(cc*cos(2*1:length(a)*pi*t/T)+d*sin(2*1:length(a)*pi*t/T))
    return(cbind(x,y))
}
