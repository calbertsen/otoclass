
#' @export

bootstrap_data <- function(covariates,group,nGroup,nGroupTest=nGroup){

    ncov <- dim(covariates)[2]

    indx <- numeric(0)
    for(i in 1:nlevels(group))
        indx <- c(indx,sample(which(as.numeric(group)==i),size=nGroup[i]))

    indxTest <- numeric(0)
    for(i in 1:nlevels(group))
        indxTest <- c(indxTest,sample(which(as.numeric(group[-indx])==i),size=nGroupTest[i]))


    means <-apply(as.matrix(covariates[indx,]),2,mean)
    sds <- apply(as.matrix(covariates[indx,]),2,function(x)diff(range(x)))
    sds <- apply(as.matrix(covariates[indx,]),2,sd)

    train <- as.matrix(covariates[indx,])
    for(i in 1:dim(train)[2]){
        train[,i] <- (train[,i]-means[i])/sds[i]
    }
    groupTrain <- factor(group[indx])
    if(ncov > 1){
        test <- as.matrix(covariates[-indx,][indxTest,])
    }else{
        test <- matrix(covariates[-indx,][indxTest],ncol=1)
    }
    for(i in 1:dim(test)[2]){
        test[,i] <- (test[,i]-means[i])/sds[i]
    }

    groupTest <- group[-indx][indxTest]

    colnames(train) <- colnames(covariates)
    rownames(train) <- rownames(covariates[indx,])
    for(i in 1:dim(train)[2]){
        train[,i] <- (train[,i]-means[i])/sds[i]
    }
    groupTrain <- factor(group[indx])
    if(ncov > 1){
        test <- as.matrix(covariates[-indx,][indxTest,])
    }else{
        test <- matrix(covariates[-indx,][indxTest],ncol=1)
    }
    for(i in 1:dim(test)[2]){
        test[,i] <- (test[,i]-means[i])/sds[i]
    }

    groupTest <- group[-indx][indxTest]

    for(i in 1:dim(train)[2]){
        train[,i] <- (train[,i]-means[i])/sds[i]
    }
    groupTrain <- factor(group[indx])
    if(ncov > 1){
        test <- as.matrix(covariates[-indx,][indxTest,])
    }else{
        test <- matrix(covariates[-indx,][indxTest],ncol=1)
    }
    for(i in 1:dim(test)[2]){
        test[,i] <- (test[,i]-means[i])/sds[i]
    }

    groupTest <- group[-indx][indxTest]

    colnames(train) <- colnames(covariates)
    rownames(train) <- names(covariates[indx])
    colnames(test) <- colnames(covariates)
    rownames(test) <- names(covariates[indx][indxTest])

    return(list(train = train,
                groupTrain = groupTrain,
                test = test,
                groupTest = groupTest
                ))

}
