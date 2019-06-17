#########################################################
######### Test read_image on different polygons #########
#########################################################

suppressPackageStartupMessages(library(rgeos))
library(sp)


efds <- read.csv("data/efd.csv",sep=",")[,-1]

picDim <- c(2500,2500)


## Test black otolith on white background

for(i in 1:nrow(efds)){
    f <- tempfile(fileext=".png")

    meanEFDList <- split(unlist(efds[i,]),gsub("[1-9][0-9]*$","",colnames(efds)))
    meanEFDList$n <- 1000
    xx <- do.call(efd2coord,meanEFDList)

    png(f,width=picDim[1],height=picDim[2])
    par(bg="white", mar = c(0,0,0,0), oma = c(0,0,0,0))
    plot.new()
    par(usr = c(1,picDim[1],1,picDim[2]))
    polygon(xx, border=NA,col="black")
    dev.off()

    im <- read_image(f)

    ## is_equal(class(im), "otolith_image")
    ## is_equal(length(im),1)
    ## is_equal(class(im[[1]]), "otolith_contour")

    
    spim <- SpatialPolygons(list(Polygons(list(Polygon(im[[1]][,1:2])),1L)),1L)
    spx <- SpatialPolygons(list(Polygons(list(Polygon(xx)),1L)),1L)
    areaX <- Reduce("+",lapply(spx@polygons,function(x)x@area))
    spdiff <- gSymdifference(spim,spx)
    is_true(Reduce("+",lapply(spdiff@polygons,function(x)x@area)) / areaX < 7 / 1000 )## / prod(picDim) < 0.1 / 1000)
}


## Test white otolith on black background

for(i in 1:nrow(efds)){
    f <- tempfile(fileext=".png")

    meanEFDList <- split(unlist(efds[i,]),gsub("[1-9][0-9]*$","",colnames(efds)))
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
    is_true(Reduce("+",lapply(spdiff@polygons,function(x)x@area)) / areaX < 7 / 1000 )## / prod(picDim) < 0.1 / 1000)
}

