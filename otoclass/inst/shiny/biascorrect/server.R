
library(shiny)
library(otoclass)
library(mvtnorm)
source("helper.R")



shinyServer(function(input,output){

    #Simulate a data set
    simData <- reactive({
        nGroup <- round(input$Ns * c(input$ps,1-input$ps))
        sigma <- list(input$sd * matrix(c(1,input$rho,input$rho,1),ncol=2),
                      matrix(c(1,input$rho,input$rho,1),ncol=2))
        mu <- cbind(c(input$x1,input$y1),c(0,0))

        tr <-array(dim=c(0,2))
        gr <- numeric(0)
        for(i in 1:length(nGroup)){
            tr <- rbind(tr,rmvnorm(sum(nGroup[i])*100,mu[,i],sigma[[i]]))
            gr <- c(gr,rep(i,sum(nGroup[i])*100))
        }
        data.frame(group=factor(gr),tr)
    })
        
    bootstrap <- reactive({
        n <- input$Nb
        
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Preparing for bootstrap", value = 0)
        Mknn1 <- list()
        Mknn3 <- list()
        Mqda <- list()
        Mlda <- list()
        dat <- simData()
        nGroup <- round(input$Ns * c(input$ps,1-input$ps))
        nGroupT <- round(input$Ns * c(input$pstest,1-input$pstest))
        indxUse <- 1:2
        resuknn1 <- array(dim=c(n,nlevels(dat[,1])))
        resuknn3 <- array(dim=c(n,nlevels(dat[,1])))
        resuqda <- array(dim=c(n,nlevels(dat[,1])))
        resulda <- array(dim=c(n,nlevels(dat[,1])))

        resknn1 <- array(dim=c(n,nlevels(dat[,1])))
        resknn3 <- array(dim=c(n,nlevels(dat[,1])))
        resqda <- array(dim=c(n,nlevels(dat[,1])))
        reslda <- array(dim=c(n,nlevels(dat[,1])))

        progress$set(message = "Bootstrapping", value = 0)
        for(r in 1:n){
            progress$inc(0, detail = paste("confusion matrix", r))
            bdat <- bootstrap_data(dat[,-1][,indxUse],
                                   factor(dat[,1]),
                                   nGroup,nGroupT)
            fitknn1 <- otoclass::knn(train=bdat$train,
                                     test=bdat$test,
                                     group=bdat$groupTrain,
                                     kn=1)
            fitknn3 <- otoclass::knn(train=bdat$train,
                                     test=bdat$test,
                                     group=bdat$groupTrain,
                                     kn=3)
            fitqda <- discrim(train=bdat$train,
                              test=bdat$test,
                              group=bdat$groupTrain,
                              type="qda",
                              verbose=FALSE)
            fitlda <- discrim(train=bdat$train,
                              test=bdat$test,
                              group=bdat$groupTrain,
                              type="lda",
                              verbose=FALSE)
            Mknn1[[r]]<-apply(table(fitknn1$predicted,bdat$groupTest),2,
                              function(x)x/sum(x))
            Mknn3[[r]]<-apply(table(fitknn3$predicted,bdat$groupTest),2,
                              function(x)x/sum(x))
            Mqda[[r]]<-apply(table(fitqda$predicted,bdat$groupTest),2,
                             function(x)x/sum(x))
            Mlda[[r]]<-apply(table(fitlda$predicted,bdat$groupTest),2,
                             function(x)x/sum(x))
  
           
            progress$inc(0.5/n, detail = paste("proportions from new sample", r))
            #New sample to test
            bdat <- bootstrap_data(dat[,-1][,indxUse],factor(dat[,1]),nGroup,nGroupT)
            fitknn1 <- otoclass::knn(train=bdat$train,
                                     test=bdat$test,
                                     group=bdat$groupTrain,
                                     kn=1)
            fitknn3 <- otoclass::knn(train=bdat$train,
                                     test=bdat$test,
                                     group=bdat$groupTrain,
                                     kn=3)
            fitqda <- discrim(train=bdat$train,
                              test=bdat$test,
                              group=bdat$groupTrain,
                              type="qda",
                              verbose=FALSE)
            fitlda <- discrim(train=bdat$train,
                              test=bdat$test,
                              group=bdat$groupTrain,
                              type="lda",
                              verbose=FALSE)

            resuknn1[r,] <- table(fitknn1$predicted)/length(fitknn1$predicted)
            resuknn3[r,] <- table(fitknn3$predicted)/length(fitknn3$predicted)
            resulda[r,] <- table(fitlda$predicted)/length(fitlda$predicted)
            resuqda[r,] <- table(fitqda$predicted)/length(fitqda$predicted)
            
            if(abs(det(Mknn1[[r]]))>1e-16)
                resknn1[r,] <- M2prop(fitknn1$predicted,list(Mknn1[[r]]))
            if(abs(det(Mknn3[[r]]))>1e-16)
                resknn3[r,] <- M2prop(fitknn3$predicted,list(Mknn3[[r]]))
            if(abs(det(Mqda[[r]]))>1e-16)
                resqda[r,] <- M2prop(fitqda$predicted,list(Mqda[[r]]))
            if(abs(det(Mlda[[r]]))>1e-16)
                reslda[r,] <- M2prop(fitlda$predicted,list(Mlda[[r]]))
            progress$inc(0.5/n)
        }
                                        
        return(list(Mknn1=Mknn1,
                    Mknn3=Mknn3,
                    Mqda=Mqda,
                    Mlda=Mlda,
                    resknn1=resknn1,
                    resknn3=resknn3,
                    resqda=resqda,
                    reslda=reslda,
                    resuknn1=resuknn1,
                    resuknn3=resuknn3,
                    resuqda=resuqda,
                    resulda=resulda)) 
    })

    output$proportionplot <- renderPlot({
        ll <- bootstrap()
        layout(matrix(1:4,ncol=2))
        if(any(!is.na(ll$resknn1[,1]))){
            plot(density(ll$resknn1[,1],na.rm=TRUE,from=0,to=1),main="1-NN",col="black")
            abline(v=input$pstest)
            abline(v=mean(ll$resknn1[,1],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resknn1[,1],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="1-NN",type="n")
        }
        if(any(!is.na(ll$resknn3[,1]))){
            plot(density(ll$resknn3[,1],na.rm=TRUE,from=0,to=1),main="3-NN",col="red")
            abline(v=input$pstest)
            abline(v=mean(ll$resknn3[,1],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resknn3[,1],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="3-NN",type="n")
        }
        if(any(!is.na(ll$resqda[,1]))){
            plot(density(ll$resqda[,1],na.rm=TRUE,from=0,to=1),main="QDA",col="orange")
            abline(v=input$pstest)
            abline(v=mean(ll$resqda[,1],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resqda[,1],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="QDA",type="n")
        }
        if(any(!is.na(ll$reslda[,1]))){
            plot(density(ll$reslda[,1],na.rm=TRUE,from=0,to=1),main="LDA",col="blue")
            abline(v=input$pstest)
            abline(v=mean(ll$reslda[,1],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$reslda[,1],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="LDA",type="n")
        }
     })

    output$rawproportionplot <- renderPlot({
        ll <- bootstrap()
        layout(matrix(1:4,ncol=2))
        if(any(!is.na(ll$resuknn1[,1]))){
            plot(density(ll$resuknn1[,1],na.rm=TRUE,from=0,to=1),main="1-NN",col="black")
            abline(v=input$pstest)
            abline(v=mean(ll$resuknn1[,1],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resuknn1[,1],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="1-NN",type="n")
        }
        if(any(!is.na(ll$resuknn3[,1]))){
            plot(density(ll$resuknn3[,1],na.rm=TRUE,from=0,to=1),main="3-NN",col="red")
            abline(v=input$pstest)
            abline(v=mean(ll$resuknn3[,1],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resuknn3[,1],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="3-NN",type="n")
        }
        if(any(!is.na(ll$resuqda[,1]))){
            plot(density(ll$resuqda[,1],na.rm=TRUE,from=0,to=1),main="QDA",col="orange")
            abline(v=input$pstest)
            abline(v=mean(ll$resuqda[,1],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resuqda[,1],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="QDA",type="n")
        }
        if(any(!is.na(ll$resulda[,1]))){
            plot(density(ll$resulda[,1],na.rm=TRUE,from=0,to=1),main="LDA",col="blue")
            abline(v=input$pstest)
            abline(v=mean(ll$resulda[,1],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resulda[,1],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="LDA",type="n")
        }
    })

    
    output$summarytable <- renderTable({
        ll <- bootstrap()
        s1 <- t(summary_wna(ll$resknn1[,1]))
        tab <- cbind(Method=c("1-NN","3-NN","QDA","LDA"),
                     rbind(s1,
                           t(summary_wna(ll$resknn3[,1])),
                           t(summary_wna(ll$resqda[,1])),
                           t(summary_wna(ll$reslda[,1])))
                     )
        #colnames(tab) <- c("Method",colnames(s1))
    })

    output$rawsummarytable <- renderTable({
        ll <- bootstrap()
        s1 <- t(summary_wna(ll$resuknn1[,1]))
        tab <- cbind(Method=c("1-NN","3-NN","QDA","LDA"),
                     rbind(s1,
                           t(summary_wna(ll$resuknn3[,1])),
                           t(summary_wna(ll$resuqda[,1])),
                           t(summary_wna(ll$resulda[,1])))
                     )
        #colnames(tab) <- c("Method",colnames(s1))
    })

    output$dataplot <- renderPlot({
         dat <- simData()
         nGroup <- round(input$Ns * c(input$ps,1-input$ps))
         nGroupT <- round(input$Ns * c(input$pstest,1-input$pstest))
         indxUse <- 1:2
         bdat <- bootstrap_data(dat[,-1][,indxUse],
                                factor(dat[,1]),
                                nGroup,nGroupT)
         
         layout(matrix(1:3,ncol=3))
         plot(dat[,2:3],col=dat[,1],main = "Synthetic data set")
         plot(bdat$train, col=bdat$groupTrain, main = "Example of bootstrapped training data")
         plot(bdat$test, col=bdat$groupTest, main = "Example of bootstrapped test data")

    })


    output$confusionmat <- renderPlot({
        ll <- bootstrap()
        layout(matrix(1:4,ncol=2,byrow=TRUE))

        plot(density(unlist(lapply(ll$Mknn1,function(x)x[1,1])),from=0,to=1,na.rm=TRUE),main="Proportion of Group 1 classified as group 1",col="black")
        lines(density(unlist(lapply(ll$Mknn3,function(x)x[1,1])),from=0,to=1,na.rm=TRUE),col="red")
        lines(density(unlist(lapply(ll$Mqda,function(x)x[1,1])),from=0,to=1,na.rm=TRUE),col="orange")
        lines(density(unlist(lapply(ll$Mlda,function(x)x[1,1])),from=0,to=1,na.rm=TRUE),col="blue")

        plot(density(unlist(lapply(ll$Mknn1,function(x)x[1,2])),from=0,to=1,na.rm=TRUE),main="Proportion of Group 2 classified as group 1")
        lines(density(unlist(lapply(ll$Mknn3,function(x)x[1,2])),from=0,to=1,na.rm=TRUE),col="red")
        lines(density(unlist(lapply(ll$Mqda,function(x)x[1,2])),from=0,to=1,na.rm=TRUE),col="orange")
        lines(density(unlist(lapply(ll$Mlda,function(x)x[1,2])),from=0,to=1,na.rm=TRUE),col="blue")

        plot(density(unlist(lapply(ll$Mknn1,function(x)x[2,1])),from=0,to=1,na.rm=TRUE),main="Proportion of Group 1 classified as group 2")
        lines(density(unlist(lapply(ll$Mknn3,function(x)x[2,1])),from=0,to=1,na.rm=TRUE),col="red")
        lines(density(unlist(lapply(ll$Mqda,function(x)x[2,1])),from=0,to=1,na.rm=TRUE),col="orange")
        lines(density(unlist(lapply(ll$Mlda,function(x)x[2,1])),from=0,to=1,na.rm=TRUE),col="blue")

        plot(density(unlist(lapply(ll$Mknn1,function(x)x[2,2])),from=0,to=1),main="Proportion of Group 2 classified as group 2")
        lines(density(unlist(lapply(ll$Mknn3,function(x)x[2,2])),from=0,to=1,na.rm=TRUE),col="red")
        lines(density(unlist(lapply(ll$Mqda,function(x)x[2,2])),from=0,to=1,na.rm=TRUE),col="orange")
        lines(density(unlist(lapply(ll$Mlda,function(x)x[2,2])),from=0,to=1,na.rm=TRUE),col="blue")

    })
    
})
