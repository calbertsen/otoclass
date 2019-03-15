#### Test calculation of elliptical fourier descriptors


## Calculation of coordinates from one harmonic
replicate(100,{
    n <- 10
    a0 <- 0
    c0 <- 0
    a <- runif(1,0,2)
    b <- runif(1,0,2)
    cc <- runif(1,0,2)
    d <- runif(1,0,2)
    xx <- efd2coord(n,a,b,cc,d,a0,c0)
    yy <- cbind(a * cos(2 * pi * (0:(n-1))/(n-1)) +
                b * sin(2 * pi * (0:(n-1))/(n-1)),
                cc * cos(2 * pi * (0:(n-1))/(n-1)) +
                d * sin(2 * pi * (0:(n-1))/(n-1)))
    is_equal(xx,yy)
})

## Normalized efd must be 1,0,0,0

n <- 10000
N <- 10000
xx <- efd2coord(n,1.0,0.0,0.0,0.5,0,0)
ee <- efd(xx,N, returnAsList = TRUE)
ee$n <- n
yy <- do.call(efd2coord,ee)


xxdx <- c(0,diff(xx[,1]))
xxdy <- c(0,diff(xx[,2]))
xxdt <- sqrt(xxdx^2 + xxdy^2)
yydx <- c(0,diff(yy[,1]))
yydy <- c(0,diff(yy[,2]))
yydt <- sqrt(yydx^2 + yydy^2)

xxt <- cumsum(xxdt) / sum(xxdt)
yyt <- cumsum(yydt) / sum(yydt)

all.equal(approx(xxt,xx[,1], xout = (0:99) / 99)$y,
          approx(yyt,yy[,1], xout = (0:99) / 99)$y)

all.equal(approx(xxt,xx[,2], xout = (0:99) / 99)$y,
          approx(yyt,yy[,2], xout = (0:99) / 99)$y)

plot(xx)
lines(yy)


