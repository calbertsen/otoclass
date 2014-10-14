#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_FACTOR(group);
  DATA_MATRIX(covariates); // #obs x #par
  DATA_MATRIX(predCovariates);
  PARAMETER_MATRIX(beta); // #par x #group-1

  Type nll;

  matrix<Type> probtmp = covariates * beta; // #obs x #group-1
  matrix<Type> prob(probtmp.rows(),probtmp.cols()+1); // #obs x #group
  prob.block(0,0,probtmp.rows(),probtmp.cols()) = probtmp;
  prob.col(prob.cols()-1) = prob.col(0)*Type(0);
  prob = exp(prob.array()).matrix();
  for(int i = 0; i < prob.rows(); ++i){
    prob.row(i) = prob.row(i)/prob.row(i).sum();
  }

  for(int i = 0; i < group.size(); i++){
    nll -= log(prob(i,group(i)));
  }

  REPORT(prob);
  ADREPORT(prob);


  matrix<Type> probtmpP = predCovariates * beta; // #obs x #group-1
  matrix<Type> probPred(probtmpP.rows(),probtmpP.cols()+1); // #obs x #group
  probPred.block(0,0,probtmpP.rows(),probtmpP.cols()) = probtmpP;
  probPred.col(probPred.cols()-1) = probPred.col(0)*Type(0);
  probPred = exp(prob.array()).matrix();//exp(probPred);
  for(int i = 0; i < probPred.rows(); ++i){
    probPred.row(i) = probPred.row(i)/probPred.row(i).sum();
  }

  REPORT(probPred);
  ADREPORT(probPred);
  

  return nll;

}

/*
library(TMB)
compile("multinom.cpp")
dyn.load(dynlib("multinom"))

##Test 1 - samme model
cova <- matrix(rnorm(1000,0,1),ncol=1)
beta <- matrix(c(1,2),ncol=2)
prob <- exp(cbind(cova %*% beta,0))
prob <- prob/rowSums(prob)
group <- array(dim=c(1000,1))
group <- sapply(1:1000,function(i){sample(c("A","B","C"),1,replace=TRUE,prob=prob[i,])})

data <- list(group = factor(group),
           covariates = cova,
	   predCovariates = matrix(seq(-5,5,0.1),ncol=1))
parameters <- list(beta = beta*0)
obj <- MakeADFun(data,parameters,DLL="multinom")
opt <- nlminb(obj$par,obj$fn,obj$gr)

## Test 2
group <- factor(sample(c("A","B"),1000,replace=TRUE,prob=c(0.5,0.5)))
grouppred <- factor(sample(c("A","B"),10000,replace=TRUE,prob=c(0.5,0.5)))
means <- (as.numeric(group)-1)*2
meanspred <- (as.numeric(grouppred)-1)*2
cova <- matrix(rnorm(1000,means,1),ncol=1)
covapred <- matrix(rnorm(10000,meanspred,1),ncol=1)

data <- list(group = factor(group),
           covariates = cbind(1,cova),
	   predCovariates = cbind(1,covapred))
parameters <- list(beta = matrix(0,nrow=2,ncol=1))
obj <- MakeADFun(data,parameters,DLL="multinom")
opt <- nlminb(obj$par,obj$fn,obj$gr)



pp<-obj$report()$probPred
plot(seq(-5,5,0.1),pp[,1],type="l",ylim=c(0,1))
lines(seq(-5,5,0.1),pp[,2],lty=1,col="red")

srep <- summary(sdreport(obj))
srep <- matrix(srep[rownames(srep)=="probPred",2],ncol=2)

lines(seq(-5,5,0.1),pp[,1]+2*srep[,1],lty=2)
lines(seq(-5,5,0.1),pp[,1]-2*srep[,1],lty=2)
lines(seq(-5,5,0.1),pp[,2]+2*srep[,2],lty=2,col="red")
lines(seq(-5,5,0.1),pp[,2]-2*srep[,2],lty=2,col="red")

 */
