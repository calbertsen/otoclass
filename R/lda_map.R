
lda_map <- function(param){
    return(list(
        logSd = factor(matrix(1:nrow(param$logSd),
            ncol = ncol(param$logSd),
            nrow = nrow(param$logSd)
            )),
        trho =  factor(matrix(1:nrow(param$trho),
            ncol = ncol(param$trho),
            nrow = nrow(param$trho)
            ))
        ))
}
