#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS


#include "../inst/include/efd.hpp"

#define install Rf_install
#define findVar Rf_findVar
#define setAttrib Rf_setAttrib
#include <TMB.hpp>


using namespace density;
using Eigen::SparseMatrix;


template<class Type>
matrix<Type> efd2coord(matrix<Type> efd, int N){
  matrix<Type> coord(2,N);
  coord.setZero();

  for(int tt = 0; tt < N; ++tt){
    Type t = (Type)tt / (Type)(N-1);
    for(int i = 0; i < efd.cols(); ++i){
      coord(0,tt) += efd(0,i) * cos(2.0 * Type(i+1) * M_PI * t);
      coord(0,tt) += efd(1,i) * sin(2.0 * Type(i+1) * M_PI * t);
      coord(1,tt) += efd(2,i) * cos(2.0 * Type(i+1) * M_PI * t);
      coord(1,tt) += efd(3,i) * sin(2.0 * Type(i+1) * M_PI * t);
    }
  }
  return coord;
}



namespace Model{
  enum Types {
    Independent,
    Unstructured,
    GMRF,
    RotatingOU
  };
}


template<class Type>
class OTOLITH_t {
  
private:
  Model::Types model;
  VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > unstructured;
  VECSCALE_t<N01<Type> > independent;
  SEPARABLE_t<GMRF_t<Type>, VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > > gmrf;

public:
  OTOLITH_t();
  OTOLITH_t(vector<Type> corpar, vector<Type> sigma, SparseMatrix<Type> Q, Model::Types model_);
  Type operator()(array<Type> x);

};

template<class Type>
OTOLITH_t<Type>::OTOLITH_t(){}

template<class Type>
OTOLITH_t<Type>::OTOLITH_t(vector<Type> corpar, vector<Type> sigma, SparseMatrix<Type> Q, Model::Types model_) :
  model(model_)
{
  switch(model) {
  case Model::Independent:
    independent = VECSCALE(N01<Type>(), sigma);
    break;
  case Model::Unstructured:
    unstructured = VECSCALE(UNSTRUCTURED_CORR(corpar), sigma);
    break;
  case Model::GMRF:
    gmrf = SEPARABLE(GMRF(Q),VECSCALE(UNSTRUCTURED_CORR(corpar), sigma));
    break;
  default:
    Rf_error("Unknown model type");
    break;
  }
}

template<class Type>
Type OTOLITH_t<Type>::operator() (array<Type> x){
  Type r = 0.0;
  switch(model) {
  case Model::Independent:
    r = independent((vector<Type>)x);
    break;
  case Model::Unstructured:
    r= unstructured((vector<Type>)x);
    break;
  case Model::GMRF:
    r= gmrf(x);
    break;
  default:
    Rf_error("Unknown model type");
    break;
  }
  return r;
}


template <class Type>
Type objective_function<Type>::operator() () {
  DATA_INTEGER(model);
  DATA_ARRAY(X);
  DATA_FACTOR(G);
  DATA_SPARSE_MATRIX(Q);

  Model::Types modelType = static_cast<Model::Types>(model);
  
  DATA_INTEGER(penalty);
  DATA_VECTOR(logLambda);// DATA_UPDATE(logLambda);
  DATA_VECTOR(prior);
  DATA_ARRAY(X_pred);
  DATA_SPARSE_MATRIX(covar);
  DATA_SPARSE_MATRIX(covar_pred);

  PARAMETER_ARRAY(mu);
  PARAMETER_ARRAY(efd);
  PARAMETER_MATRIX(logSigma);
  PARAMETER_MATRIX(corpar);
  PARAMETER_VECTOR(logDelta);


  Type scale = X.dim.prod();
  
  vector<Type> lambda = logLambda.exp();
  vector<Type> delta = logDelta.exp();

  matrix<Type> sigma(logSigma.rows(), logSigma.cols());
  for(int i = 0; i < sigma.rows(); ++i)
    for(int j = 0; j < sigma.cols(); ++j)
      sigma(i,j) = exp(logSigma(i,j));

  vector<SparseMatrix<Type> > Qg(NLEVELS(G));
  for(int j = 0; j < NLEVELS(G); ++j){
    Qg(j) = Q;
    for(int i = 0; i < Q.cols(); ++i){
      Qg(j).coeffRef(i,i) += delta(j);
    }
  }
  REPORT(Qg(0));

  vector<matrix<Type> > coords(NLEVELS(G));
  matrix<Type> meanCoords(2,X.dim[1]);
  meanCoords.setZero();
  if(efd.cols() > 0){
    for(int i = 0; i < NLEVELS(G); ++i){
      coords(i) = efd2coord(efd.col(i).matrix(), X.dim[1]); // Add something to account for length etc.
      meanCoords += coords(i);
    }
    meanCoords /= (Type)NLEVELS(G);
  }

  // Estimate
  vector<OTOLITH_t<Type> > dist(NLEVELS(G));
  for(int i = 0; i < NLEVELS(G); ++i)
    dist(i) = OTOLITH_t<Type>((vector<Type>)corpar.col(i), (vector<Type>)sigma.col(i), Qg(i), modelType);
  
  Type nll = 0.0;

  for(int i = 0; i < X.cols(); ++i){
    array<Type> tmp = X.col(i);
    switch(modelType) {
    case Model::Independent:
      tmp -= (vector<Type>)(covar.col(i).transpose() * mu.col(G(i)).matrix());
      break;
    case Model::Unstructured:
      tmp -= (vector<Type>)(covar.col(i).transpose() * mu.col(G(i)).matrix());
      break;
    case Model::GMRF:
      tmp -= coords(G(i)).array();
      break;
    default:
      Rf_error("Unknown model type");
      break;
    }

    nll += dist(G(i))(tmp) / scale;
  }


  // Penalty
  // vector<Type> muRMean = mu.rowwise().mean();
  if(penalty > 0){

    array<Type> muCMean(mu.dim[0],mu.dim[1]); // Mean over groups
    muCMean.setZero();
    for(int k = 0; k < mu.dim[2]; ++k)
      muCMean += mu.col(k);
    muCMean /= mu.dim[2];

    REPORT(muCMean);
    
    Type normR = 0.0;
    Type normC = 0.0;

    //Type normC = 0.0;
    // Penalty if mu is used
    for(int i = 0; i < mu.dim[0]; ++i)
      for(int j = 0; j < mu.dim[1]; ++j){
	normC += pow(abs(muCMean(i,j)), penalty);
	for(int k = 0; k < mu.dim[2]; ++k){
	  normR += pow(abs(mu(i,j,k) - muCMean(i,j)), penalty); // sigma(i,j) *
	  //normR += pow(abs(mu(i,j,k)), penalty); // sigma(i,j) *
	}
      }
    nll += lambda(0) * pow(normR, 1.0/(double)penalty);
    nll += lambda(1) * pow(normC, 1.0/(double)penalty);

    Type normCor = 0.0;
    for(int j = 0; j < corpar.cols(); ++j)
      for(int i = 0; i < corpar.rows(); ++i)
    	normCor += pow(abs(corpar(i,j)), penalty);
    nll += lambda(2) * pow(normCor, 1.0/(double)penalty);
    
    // Penalty if efd is used
    if(efd.cols() > 0){
      for(int j = 0; j < NLEVELS(G); ++j){
	// for(int cc = 0; cc < meanCoords.rows(); ++cc){
	//   Type norm = 0.0;
	//   for(int i = 0; i < meanCoords.cols(); ++i){
	//     norm += pow(abs(coords(j)(cc,i) - meanCoords(cc,i)), penalty); 
	//   }
	//   nll += lambda * pow(norm, 1.0/(double)penalty);	  
	// }
	for(int cc = 0; cc < efd.dim[0]; ++cc){
	  Type norm = 0.0;
	  for(int i = 1; i < efd.dim[1]; ++i){
	    norm += pow(abs(efd(cc,i,j)), penalty); 
	  }
	  if(cc == 3)
	    norm += pow(abs(efd(cc,0,j)), penalty); 
	  nll += lambda(0) * pow(norm, 1.0/(double)penalty);	  
	}
      }
    }
  }

  // Predict
  matrix<Type> logpred(NLEVELS(G),X_pred.cols());
  for(int i = 0; i < logpred.cols(); ++i){
    vector<Type> lp(NLEVELS(G));
    Type lps = 0.0;
    for(int j = 0; j < NLEVELS(G); ++j){
      array<Type> tmpPred;
      tmpPred = X_pred.col(i);
      switch(modelType) {
      case Model::Independent:
        tmpPred -= (vector<Type>)(covar_pred.col(i).transpose() * mu.col(j).matrix());
	break;
      case Model::Unstructured:
	tmpPred -= (vector<Type>)(covar_pred.col(i).transpose() * mu.col(j).matrix());
	break;
      case Model::GMRF:
        tmpPred -= coords(j).array();
	break;
      default:
	Rf_error("Unknown model type");
	break;
      }
      lp(j) = -dist(j)(tmpPred) + log(prior(j));
      if(j == 0){
	lps = lp(j);
      }else{
	lps = logspace_add(lps,lp(j));
      }
    }
    logpred.col(i) = lp - lps;
  }
  REPORT(logpred);
  REPORT(coords);
  REPORT(meanCoords);
  
  return nll;
}
