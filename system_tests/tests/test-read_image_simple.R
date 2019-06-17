#########################################################
######### Test read_image on different polygons #########
#########################################################

suppressPackageStartupMessages(library(rgeos))
library(sp)

set.seed(1)

meanEFD <- c(A1 = 435.735280408116, A2 = 1.23083957873596,
             B1 = 5.96353787417546, B2 = 4.90385332405064,
             C1 = -9.63023759177977, C2 = -23.6630393747043,
             D1 = 231.227847131746, D2 = -2.5440249347623,
             A0 = 815.481837821778,  
             C0 = 550.369870693356)

corpar <- c(0.156315014270222, -0.000693091825817576, 0.030545388214153, 0.209692848472478,
            0.156248861793837, -0.692155670981122, -0.141311060323548, 0.18322919679499,
            -0.732753181778909, 0.581885353974751, -1.29232001798667, -0.152731252073664,
            0.567329068085384, -0.795209180219285, -0.0730509111752459, 3.04111933832025,
            0.0392478334348732, -0.0458353646295011, 0.0742861398235585, -0.217660965285963,
            -0.449401507986472, -0.0438991342766151, 0.151434197999488, -0.551890654377699,
            0.163442762464386, 0.306931450016406, -0.169347566553431, -0.0827722863765649,
            0.0707633335815199, -0.0554823403104922, -0.086566127378944, -0.0547408043852815,
            -0.081963816620411, 0.0334181051629354, 0.0064502488425685, 0.024031196458353,
            0.0264536148215254, -0.0173111088621229, 0.000109104942350734, -0.0117190629870908,
            0.121542950520781, -0.0722354280102189, -0.0988993401038516, 0.0531498802869901,
            0.393048684142256)

sigma <- diag(exp(c(4.08074626456868, 1.67711538307622, 3.31370811676143, 2.49041812851168,
           2.74542957227797, 2.46750394555354, 3.69596942268887, 1.86307575436782,
           0.1*4.33065081109118, 0.1*5.8304339151742)))

Cm <- diag(1,length(meanEFD))
Cm[lower.tri(Cm)] <- corpar
Sigma <- sigma %*% cov2cor(Cm %*% t(Cm)) %*% t(sigma)


simEFD <- function(){
    xx <- as.vector(meanEFD + t(chol(Sigma)) %*% rnorm(length(meanEFD)))
    names(xx) <- names(meanEFD)
    return(xx)
    }


picDim <- c(1500,1500)


## Test black otolith on white background

invisible(
    replicate(100,{
        f <- tempfile(fileext=".png")

        meanEFDList <- split(simEFD(),gsub("[1-9][0-9]*$","",names(meanEFD)))
        meanEFDList$n <- 1000
        xx <- do.call(efd2coord,meanEFDList)

        png(f,width=picDim[1],height=picDim[2])
        par(bg="white", mar = c(0,0,0,0), oma = c(0,0,0,0))
        plot.new()
        par(usr = c(1,picDim[1],1,picDim[2]))
        polygon(xx, border=NA,col="black")
        dev.off()

        im <- read_image(f)
        is_equal(class(im), "otolith_image")
        is_equal(length(im),1)
        is_equal(class(im[[1]]), "otolith_contour")

        spim <- SpatialPolygons(list(Polygons(list(Polygon(im[[1]][,1:2])),1L)),1L)
        spx <- SpatialPolygons(list(Polygons(list(Polygon(xx)),1L)),1L)
        areaX <- Reduce("+",lapply(spx@polygons,function(x)x@area))
        spdiff <- gSymdifference(spim,spx)
        is_true(Reduce("+",lapply(spdiff@polygons,function(x)x@area)) / areaX < 6 / 1000 )##prod(picDim) < 0.1 / 1000)
    })
)


## Test white otolith on black background

invisible(
    replicate(100,{
        f <- tempfile(fileext=".png")

        meanEFDList <- split(simEFD(),gsub("[1-9][0-9]*$","",names(meanEFD)))
        meanEFDList$n <- 1000
        xx <- do.call(efd2coord,meanEFDList)

        png(f,width=picDim[1],height=picDim[2])
        par(bg="black", mar = c(0,0,0,0), oma = c(0,0,0,0))
        plot.new()
        par(usr = c(1,picDim[1],1,picDim[2]))
        polygon(xx, border=NA,col="white")
        dev.off()

        im <- read_image(f)
        is_equal(class(im), "otolith_image")
        is_equal(length(im),1)
        is_equal(class(im[[1]]), "otolith_contour")

        spim <- SpatialPolygons(list(Polygons(list(Polygon(im[[1]][,1:2])),1L)),1L)
        spx <- SpatialPolygons(list(Polygons(list(Polygon(xx)),1L)),1L)
        areaX <- Reduce("+",lapply(spx@polygons,function(x)x@area))
        spdiff <- gSymdifference(spim,spx)
        is_true(Reduce("+",lapply(spdiff@polygons,function(x)x@area)) / areaX < 6 / 1000 )## / prod(picDim) < 0.1 / 1000)
    })
)

## Illustrate allowed difference
if(FALSE){
    ##png("allowed_diff.png",1000,1000)
    par(bg="white", mar = c(0,0,0,0), oma = c(0,0,0,0))
    plot.new()
    par(usr = c(1,100,1,100))
    polygon(c(0,sqrt(6),sqrt(6),0,0) + 50-0.5,
            c(0,0,sqrt(6),sqrt(6),0) + 50-0.5,
            border=NA,col="black")
    box()
    ##dev.off()
}

