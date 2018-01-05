##' Run shiny example included in the package
##'
##' @param name Name of the shiny example (biascorrect, or proportion)
##' @return Nothing
##' @author Christoffer Moesgaard Albertsen
runShinyExample <- function(name){
    if(requireNamespace("shiny")){
        homedir <- system.file("shiny", package = "otoclass")
        shiny::runApp(paste(homedir,name,sep="/"))
    }else{
        warning("The shiny package is needed to run the example")
    }
}

