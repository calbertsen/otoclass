####################################################################
######### Compare EFD calculations to the R package shapeR #########
######### https://github.com/lisalibungan/shapeR           #########
####################################################################


## Load packages
library(shapeR)

set.seed(1)

## Test EFD calculation
invisible(
    replicate(400,{
            x <- efd2coord(20,
                           runif(50,0,2),
                           runif(50,0,2),
                           runif(50,0,2),
                           runif(50,0,2),
                           0,0)
            e1l <- shapeR:::.shapeR.efourier(x, n = 20L)
            e1l$ao <- e1l$co <- NULL
            e1 <- do.call(cbind,e1l)
            e2l <- efd(x, 20, returnAsList = TRUE)
            e2l$A0 <- e2l$C0 <- NULL
            e2 <- do.call(cbind,e2l)
            is_equal(e1,e2)
            e3l <- shapeR:::.shapeR.NEF(x, n = 20L)
            e3l$ao <- e3l$co <- e3l$size <- e3l$theta <- e3l$psi <- NULL
            e3 <- do.call(cbind,e3l)
            e4l <- efd(x, 20, returnAsList = TRUE, normalize = TRUE)
            e4l$A0 <- e4l$C0 <- NULL
            e4 <- do.call(cbind,e4l)
            is_equal(e3,e4)
    })
)


## Test efd2coord

invisible(
    replicate(400,{
            x <- efd2coord(20,
                           runif(50,0,2),
                           runif(50,0,2),
                           runif(50,0,2),
                           runif(50,0,2),
                           0,0)
            ee <- efd(x,20, returnAsList = TRUE)
            ee$n <- 200
            ee2 <- ee
            ee2$n <- ee2$n - 1
            ee2$A0 <- ee2$A0 * 2
            ee2$C0 <- ee2$C0 * 2
            names(ee2) <- c("an","bn","cn","dn","ao","co","n")
            ee2$k <- length(ee2$an)
            xx1l <- do.call(shapeR:::.shapeR.iefourier,ee2)
            xx1 <- do.call(cbind,xx1l)
            xx2 <- do.call(efd2coord,ee)
            is_equal(xx1,xx2[-nrow(xx2),])
    })
)
