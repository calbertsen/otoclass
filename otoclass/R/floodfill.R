
floodfill <- function(pic, tolerance = 0, x0 = 0, y0 = 0, setInitCol = FALSE){
    pic <- round(.Call("scanlineFill", pic/255, tolerance, as.integer(x0), as.integer(y0), setInitCol) * 
                 255)
    pic[pic == -255] <- 0
    pic
}
