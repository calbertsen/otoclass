#######################################################################
######### Compare EFD calculations to the Python module pyefd #########
######### https://github.com/hbldh/pyefd                      #########
#######################################################################


## Load packages
library(reticulate)

set.seed(1)

## Test EFD calculation
invisible(
    replicate(400,{
            numpy <- import("numpy")
            pyefd <- import("pyefd")
            x <- efd2coord(20,
                           runif(50,0,2),
                           runif(50,0,2),
                           runif(50,0,2),
                           runif(50,0,2),
                           0,0)
            e1 <- pyefd$elliptic_fourier_descriptors(x, order = 20L)
            e2l <- efd(x, 20, returnAsList = TRUE)
            e2l$A0 <- e2l$C0 <- NULL
            e2 <- do.call(cbind,e2l)
            is_equal(e1,e2)
            e3 <- pyefd$elliptic_fourier_descriptors(x, order = 20L, normalize = TRUE)
            e4l <- efd(x, 20, returnAsList = TRUE, normalize = TRUE)
            e4l$A0 <- e4l$C0 <- NULL
            e4 <- do.call(cbind,e4l)
            is_equal(e3,e4)
    })
)
