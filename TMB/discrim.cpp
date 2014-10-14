#include <TMB.hpp>

using namespace density;

template <class Type>
Type objective_function<Type>::operator()()
{
  DATA_FACTOR(group);
  DATA_MATRIX(covariates); // #obs x #par
  DATA_VECTOR(pi); // #group : prior prob of belonging to group
  DATA_INTEGER(distribution); // Distribution of covariates given group
  DATA_MATRIX(predCovariates);
  PARAMETER_MATRIX(means); // #par x #group
  PARAMETER_MATRIX(logSd); // #par x #group
  PARAMETER_MATRIX(trho); // #par*(#par-1)/2 x #group

  Type nll = 0.0;

  array<Type> covs(covariates.cols(),covariates.cols(),means.cols()); // #par x # par # # group

  int nrho = 0;
  for(int k = 0; k < means.cols(); ++k){ // Grupper
    for(int i = 0; i < covariates.cols(); ++i){ // Par
      covs(i,i,k) = exp(2.0*logSd(i,k));
      for(int j = i + 1; j < covariates.cols(); ++j){ //Resten af par
	Type rho = trho(nrho,k)/sqrt(Type(1.0)+trho(nrho,k)*trho(nrho,k));
	covs(i,j,k) = rho*exp(logSd(i,k)+logSd(j,k));
	covs(j,i,k) = covs(i,j,k);
	nrho += 1;
      }
    }
    nrho = 0;
  }
  

  //Train
  matrix<Type> prob(covariates.rows(),means.cols());

  for(int i = 0; i < prob.rows(); ++i){
    for(int k = 0; k < means.cols(); ++k){
      vector<Type> tmp(means.rows());
      tmp = covariates.row(i);
      tmp -= (vector<Type>)means.col(k);
      array<Type> tmpcovmat = covs.col(k);
      matrix<Type> covmat(tmpcovmat.dim(0),tmpcovmat.dim(1));
      for(int indi = 0; indi < covmat.rows(); ++indi){
	for(int indj = 0; indj < covmat.cols(); ++indj){
	  covmat(indi,indj) = tmpcovmat(indi,indj);
	}
      }
      prob(i,k) = exp(-MVNORM(covmat)(tmp));
      if(k == group(i)){
	nll += MVNORM(covmat)(tmp);
      }
    }
    prob.row(i) = prob.row(i)/prob.row(i).sum();

    //nll -= log(prob(i,group(i)));
  }

  REPORT(prob);
  //ADREPORT(prob);
  
  
  //Test
  matrix<Type> probPred(predCovariates.rows(),means.cols());

  for(int i = 0; i < probPred.rows(); ++i){
    for(int k = 0; k < means.cols(); ++k){
      vector<Type> tmp(means.rows());
      tmp = predCovariates.row(i);
      tmp -= (vector<Type>)means.col(k);
      array<Type> tmpcovmat = covs.col(k);
      matrix<Type> covmat(tmpcovmat.dim(0),tmpcovmat.dim(1));
      for(int indi = 0; indi < covmat.rows(); ++indi){
	for(int indj = 0; indj < covmat.cols(); ++indj){
	  covmat(indi,indj) = tmpcovmat(indi,indj);
	}
      }
      probPred(i,k) = exp(-MVNORM(covmat)(tmp));     
    }
    probPred.row(i) = probPred.row(i)/probPred.row(i).sum();
  }

  REPORT(probPred);
  //ADREPORT(probPred);
  
  REPORT(covs);

  

  return nll;

}
