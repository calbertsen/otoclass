
library(shiny)
library(otoclass)
library(mvtnorm)
source('helper.R')

dattodown <- read.csv("example.csv")


Sys.setlocale(locale="C")
options(shiny.maxRequestSize=60*1024^2)

shinyServer(function(input,output){

    output$downloadData <- downloadHandler(
        filename = "example.csv",
        content = function(file) {
            write.csv(dattodown,file, row.names=FALSE)
        }
    )
    
    
    readDat <- reactive({
        inFile <- input$datfile
        if (is.null(inFile))
            return(NULL)
        read.csv(inFile$datapath, header = input$header,
                 sep = input$sep, quote = input$quote)
    })
    output$selectShowGroup <- renderUI({
        dat <- readDat()
        if (is.null(dat))
            return(NULL)
        group <- factor(dat[,input$classgroup])
        choice <- 1:nlevels(group)
        names(choice) <- levels(group)
        
        selectInput("showgrp","Select group to show:", choices = choice) 
    })
     output$selectVariables <- renderUI({
        dat <- readDat()
        if (is.null(dat))
            return(NULL)
        nam <- colnames(dat)
        #choice <- 1:length(nam)
        #names(choice) <- nam
        if(!is.null(dat))
        checkboxGroupInput("indxUse", "Choose the class label variable:", nam,inline=TRUE)
    })
        output$selectGroup <- renderUI({
        dat <- readDat()
        if (is.null(dat))
            return(NULL)
        nam <- colnames(dat)
        #choice <- 1:length(nam)
        #names(choice) <- nam
        if(!is.null(dat))
        selectInput("classgroup","Choose the predictors:", nam)
    })

    output$slidersTrain <- renderUI({
        dat <- readDat()
        if (is.null(dat))
            return(NULL)
        group <- factor(dat[,input$classgroup])
        numprgrp <- table(group)
        namegrp <- levels(group)
        numslide <- nlevels(group)
        lapply(1:numslide, function(i) {
            numericInput(inputId = paste0("ngroupTrain", i), label = paste0(namegrp[i]," (max: ",floor(numprgrp[i]*0.5),"):"),
                        min = 1, max = floor(numprgrp[i]*0.5), value = min(floor(numprgrp*0.5)), step = 1)
        })
    })
    output$slidersTest <- renderUI({
        dat <- readDat()
        if (is.null(dat))
            return(NULL)
        group <- factor(dat[,input$classgroup])
        numprgrp <- table(group)
        namegrp <- levels(group)
        numslide <- nlevels(group)
        lapply(1:numslide, function(i) {
           numericInput(inputId = paste0("ngroupTest", i), label = paste0(namegrp[i]," (max: ",floor(numprgrp[i]*0.5),"): "),
                        min = 1, max = floor(numprgrp[i]*0.5), value = min(floor(numprgrp*0.5)), step = 1)
 
        })

    })

    output$uploadHead <- renderTable({
        head(readDat())
    })
    output$groupTable <- renderTable({
        if (is.null(readDat()))
            return(NULL)
 
        t(as.matrix(table(factor(readDat()[,input$classgroup]))))
    })
    output$plotVars <- renderPlot({
        dat <- readDat()
        if (is.null(dat))
            return(NULL)
        if(is.null(input$indxUse))
            return(NULL)
        group <- factor(dat[,input$classgroup])
        if(length(input$indxUse)>1){
            if(length(input$indxUse)>6){
                pairs(dat[,input$indxUse[1:6]],col=group,main="Only the 6 first variables are shown")
            }else{
                pairs(dat[,input$indxUse],col=group,main="Only the 10 first variables are shown")
            }
        }else{
            plot(dat[,input$indxUse],group,col=group)
        }

    })

    getProp <- reactive({

        dat <- readDat()
        if(is.null(dat))
            return(NULL)
        group <- factor(dat[,input$classgroup])
        
        train <- sapply(1:nlevels(group), function(i) {
            as.numeric(input[[paste0("ngroupTrain", i)]])
        })
        test <- sapply(1:nlevels(group), function(i) {
            as.numeric(input[[paste0("ngroupTest", i)]])
        })

        ngrp <- data.frame(train=unlist(train),test=unlist(test))
        rownames(ngrp) <- levels(group)

        list(numbers = ngrp,
             proportions = apply(ngrp,2,function(x)x/sum(x))
             )
    })
         

    output$bootprop <- renderTable({
        gp <- getProp()
        gp$proportions
    })
           
    bootstrap <- reactive({
        n <- input$Nb

        nboot <- getProp()
        if(is.null(nboot))
            return(NULL)
        if(is.null(input$indxUse))
            return(NULL)
        nGroup <- nboot$numbers[,1]
        nGroupT <- nboot$numbers[,2]
        proptest <- nboot$proportions[,2]
       
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        progress$set(message = "Preparing for bootstrap", value = 0)

        dat <- readDat()
        group <- factor(dat[,input$classgroup])
        datrest <- dat[,input$indxUse]
        npred <- length(input$indxUse)	
	if(npred == 1){
           datrest <- matrix(datrest,ncol=1)
        }
        
        resknn1 <- array(dim=c(n,nlevels(group)))
        resknn3 <- array(dim=c(n,nlevels(group)))
        resqda <- array(dim=c(n,nlevels(group)))
        reslda <- array(dim=c(n,nlevels(group)))
        correctClass <- array(dim=c(n,4))
        

        progress$set(message = "Bootstrapping", value = 0)
        for(r in 1:n){
            progress$inc(0, detail = r)
            bdat <- bootstrap_data(datrest,
                                   group,
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
            
            resknn1[r,] <- table(fitknn1$predicted)/length(bdat$groupTest)
            resknn3[r,] <- table(fitknn3$predicted)/length(bdat$groupTest)
            resqda[r,] <- table(fitqda$predicted)/length(bdat$groupTest)
            reslda[r,] <- table(fitlda$predicted)/length(bdat$groupTest)

            correctClass[r,1] <- sum(diag(table(fitknn1$predicted,bdat$groupTest)))/length(bdat$groupTest)
            correctClass[r,2] <- sum(diag(table(fitknn3$predicted,bdat$groupTest)))/length(bdat$groupTest)
            correctClass[r,3] <- sum(diag(table(fitqda$predicted,bdat$groupTest)))/length(bdat$groupTest)
            correctClass[r,4] <- sum(diag(table(fitlda$predicted,bdat$groupTest)))/length(bdat$groupTest)


            progress$inc(1/n)
        }
                                        
        return(list(resknn1=resknn1,
                    resknn3=resknn3,
                    resqda=resqda,
                    reslda=reslda,
                    correctClass=correctClass)) 
    })

    output$proportionplot <- renderPlot({
        ll <- bootstrap()
        if(is.null(ll))
            return(NULL)
        nboot <- getProp()
        proptest <- nboot$proportions[,2]
        dat <- readDat()
        group <- factor(dat[,input$classgroup])
        
        indx <- as.numeric(input$showgrp)
        
        layout(matrix(1:4,ncol=2))
        if(any(!is.na(ll$resknn1[,indx]))){
            plot(density(ll$resknn1[,indx],na.rm=TRUE,from=0,to=1),main="1-NN",col="black")
            abline(v=proptest[indx])
            abline(v=mean(ll$resknn1[,indx],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resknn1[,indx],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="1-NN",type="n")
        }
        if(any(!is.na(ll$resknn3[,indx]))){
            plot(density(ll$resknn3[,indx],na.rm=TRUE,from=0,to=1),main="3-NN",col="red")
            abline(v=proptest[indx])
            abline(v=mean(ll$resknn3[,indx],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resknn3[,indx],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="3-NN",type="n")
        }
        if(any(!is.na(ll$resqda[,indx]))){
            plot(density(ll$resqda[,indx],na.rm=TRUE,from=0,to=1),main="QDA",col="orange")
            abline(v=proptest[indx])
            abline(v=mean(ll$resqda[,indx],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$resqda[,indx],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="QDA",type="n")
        }
        if(any(!is.na(ll$reslda[,1]))){
            plot(density(ll$reslda[,indx],na.rm=TRUE,from=0,to=1),main="LDA",col="blue")
            abline(v=proptest[indx])
            abline(v=mean(ll$reslda[,indx],na.rm=TRUE),lty=2,col="grey")
            abline(v=median(ll$reslda[,indx],na.rm=TRUE),lty=3,col="grey")
        }else{
            plot(NA,ylim=c(0,1),xlim=c(0,1),main="LDA",type="n")
        }
     })

 
    
    output$summarytable <- renderTable({
        ll <- bootstrap()
        if(is.null(ll))
            return(NULL)
        dat <- readDat()
        group <- factor(dat[,input$classgroup])
        
        indx <- as.numeric(input$showgrp)
        
        s1 <- t(summary_wna(ll$resknn1[,indx]))
        tab <- cbind(Method=c("1-NN","3-NN","QDA","LDA"),
                     rbind(s1,
                           t(summary_wna(ll$resknn3[,indx])),
                           t(summary_wna(ll$resqda[,indx])),
                           t(summary_wna(ll$reslda[,indx])))
                     )
        #colnames(tab) <- c("Method",colnames(s1))
    })

    output$ccplot2tal <- renderPlot({
        ll <- bootstrap()
        if(is.null(ll))
            return(NULL)
        correctClass <- ll$correctClass
        knn1 <- density(correctClass[,1],na.rm=TRUE,from=0,to=1)
        knn3 <- density(correctClass[,2],na.rm=TRUE,from=0,to=1)
        qda <- density(correctClass[,3],na.rm=TRUE,from=0,to=1)
        lda <- density(correctClass[,4],na.rm=TRUE,from=0,to=1)
        yrng <- c(0,max(c(knn1$y,knn3$y,qda$y,lda$y)))
        plot(knn1,main="Classification success (total)",col="black",ylim=yrng)
        lines(knn3,col="red")
        lines(qda,col="orange")
        lines(lda,col="blue")
    })
})
