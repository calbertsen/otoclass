#' @useDynLib otoclass
#' @importFrom Rcpp evalCpp

.onLoad <- function(libname,pkgname){
    #path <- paste(libname,pkgname,"TMB","",sep="/")
    #if(!file.exists(TMB::dynlib(paste0(path,"multinom"))))
    #    TMB::compile(paste0(path,"multinom.cpp"))
    #    if(!file.exists(TMB::dynlib(paste0(path,"discrim"))))
    #    TMB::compile(paste0(path,"discrim.cpp"))

    #dyn.load(TMB::dynlib(paste0(path,"multinom")))
    #dyn.load(TMB::dynlib(paste0(path,"discrim")))

}
