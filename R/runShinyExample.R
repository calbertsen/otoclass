
#' @export


runShinyExample <- function(name){
    require(shiny)
    homedir <- system.file("shiny", package = "otoclass")
    runApp(paste(homedir,name,sep="/"))
}

