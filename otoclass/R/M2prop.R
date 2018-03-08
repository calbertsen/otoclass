#' @export

M2prop <- function(pred,Mlist){

    ml <- lapply(Mlist,function(M)try(solve(M,table(pred)/length(pred))))
    me <- matrix(unlist(ml),ncol=length(table(pred)),byrow=TRUE)
    indxUse<-c()
    for(i in 1:dim(me)[1]){
        if(any(me[i,] < 0)){
            while(any(me[i,] < 0)){
                mis <- which(me[i,] <= 0)
                if(length(me[i,])-length(mis)>1){
                    newprop <- me[i,]
                    newprop[mis] <- 0
                    mnew <- Mlist[[i]][,-mis]
                    #mnew <- M[,-mis]
                    muse <- solve(t(mnew)%*%mnew)%*%t(mnew)
                    
                    newprop[-mis] <- muse%*%(table(pred)/length(pred))
                        
                }else{
                    newprop <- 0*me[i,]
                    newprop[-mis] <- 1
                }
                me[i,] <- newprop
            }
        }
    }
    me
}

