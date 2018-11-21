#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS


#include "../inst/include/efd.hpp"

#define install Rf_install
#define findVar Rf_findVar
#define setAttrib Rf_setAttrib
#include <TMB.hpp>





namespace my_atomic {
   template<class Float>
  struct Lp_t {
    typedef Float Scalar;
    Float lambda;
    Float p;
    Float xx;
    Float operator() (Float x) {
      return exp(-pow(fabs(x),p) / lambda) / lambda;
    }
    Float density() {
      using gauss_kronrod::integrate;
      Float norm = integrate(*this, -INFINITY, INFINITY);
      Float ans = -pow(fabs(xx),p) / lambda - log(lambda) - log(norm);
      return ans;
    }
  };

  template<class Float>
  Float eval(Float x, Float p, Float lambda) {
    Lp_t<Float> f = {x, p, lambda};
    return f.density();
  }
  TMB_BIND_ATOMIC(func, 101, eval(x[0], x[1], x[2]))
  template<class Type>
  Type Lp(Type x, Type p, Type lambda) {
    vector<Type> args(4); // Last index reserved for derivative order
    args << x, p, lambda, 0;
    return my_atomic::func(CppAD::vector<Type>(args))[0];
  }
}









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
  DATA_INTEGER(estimateUnbiasedTestProportions);
  DATA_INTEGER(trainContributionToTestProportions);
  DATA_FACTOR(testBiasCorrectionGroup);
  DATA_FACTOR(trainBiasCorrectionGroup);
  DATA_INTEGER(addMisfitCategory);
  
  Model::Types modelType = static_cast<Model::Types>(model);
  
  DATA_INTEGER(penalty);
  DATA_VECTOR(prior);
  DATA_ARRAY(X_pred);
  DATA_SPARSE_MATRIX(covar);
  DATA_SPARSE_MATRIX(covar_pred);
  DATA_SPARSE_MATRIX(commonCovar);
  DATA_SPARSE_MATRIX(commonCovar_pred);

  
  PARAMETER_ARRAY(mu);
  PARAMETER_MATRIX(commonMu);
  PARAMETER_ARRAY(efd);
  PARAMETER_MATRIX(logSigma);
  PARAMETER_MATRIX(corpar);
  PARAMETER_VECTOR(logDelta);
  PARAMETER_VECTOR(logLambda);// DATA_UPDATE(logLambda);
  PARAMETER_MATRIX(thetaIn);

  Type scale = 1.0; //X.dim.prod();
  
  vector<Type> lambda = exp(-logLambda);
  vector<Type> delta = logDelta.exp();


  array<Type> muUse(mu.dim);
    for(int i = 0; i < muUse.dim[0]; ++i)
      for(int j = 0; j < muUse.dim[1]; ++j){
	muUse(i,j,0) = mu(i,j,0);
	for(int k = 1; k < muUse.dim[2]; ++k){
	  muUse(i,j,k) = muUse(i,j,0) + mu(i,j,k);
	}
      }
  
  
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
      tmp -= (vector<Type>)(covar.col(i).transpose() * muUse.col(G(i)).matrix());
      if(commonCovar.rows() > 0){
	tmp -= (vector<Type>)(commonCovar.col(i).transpose() * commonMu);
      }
      break;
    case Model::Unstructured:
      tmp -= (vector<Type>)(covar.col(i).transpose() * muUse.col(G(i)).matrix());
      if(commonCovar.rows() > 0){
	tmp -= (vector<Type>)(commonCovar.col(i).transpose() * commonMu);
      }
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
  if(penalty != 0){

    // array<Type> muCMean(mu.dim[0],mu.dim[1]); // Mean over groups
    // muCMean.setZero();
    // for(int k = 0; k < mu.dim[2]; ++k)
    //   muCMean += muUse.col(k);
    // muCMean /= mu.dim[2];

    // REPORT(muCMean);
    
    // Type normR = 0.0;
    // Type normC = 0.0;

    //Type normC = 0.0;
    // Penalty if mu is used
    for(int i = 0; i < mu.dim[0]; ++i)
      for(int j = 0; j < mu.dim[1]; ++j){
	// nll += pow(abs(muCMean(i,j)), penalty) * lambda(1);
	for(int k = 0; k < mu.dim[2]; ++k){
	  // nll += pow(abs(mu(i,j,k) - muCMean(i,j)), penalty) * lambda(0); // sigma(i,j) *
	  switch(penalty){
	  case -1:
	    nll -= dt(mu(i,j,k) / lambda(0),Type(3.0),true) - log(lambda(0));
	    break;
	  case 1:
	    nll -= -fabs(mu(i,j,k)) / lambda(0) - log(2.0 * lambda(0));
	    break;
	  case 2:
	    nll -= dnorm(mu(i,j,k),Type(0.0),lambda(0), true); // sigma(i,j) *
	    break;
	  default:
	    nll -= my_atomic::Lp(mu(i,j,k),(Type)penalty,lambda(0)); //pow(abs(mu(i,j,k)), penalty) / lambda(0) - log(lambda(0));
	  }
	}
      }
    for(int i = 0; i < commonMu.rows(); ++i)
      for(int j = 0; j < commonMu.cols(); ++j){
	// nll += pow(abs(muCMean(i,j)), penalty) * lambda(1);
	// nll += pow(abs(mu(i,j,k) - muCMean(i,j)), penalty) * lambda(0); // sigma(i,j) *
	switch(penalty){
	case -1:
	  nll -= dt(commonMu(i,j) / lambda(0),Type(3.0),true) - log(lambda(0));
	  break;
	case 1:
	  nll -= -fabs(commonMu(i,j)) / lambda(0) - log(2.0 * lambda(0));
	  break;
	case 2:
	  nll -= dnorm(commonMu(i,j),Type(0.0),lambda(0), true); // sigma(i,j) *
	  break;
	default:
	  nll -= my_atomic::Lp(commonMu(i,j),(Type)penalty,lambda(0)); //pow(abs(mu(i,j,k)), penalty) / lambda(0) - log(lambda(0));
	}
      }
    // nll += lambda(0) * pow(normR, 1.0/(double)penalty);
    // nll += lambda(1) * pow(normC, 1.0/(double)penalty);

    // // Type normCor = 0.0;
    // for(int j = 0; j < corpar.cols(); ++j)
    //   for(int i = 0; i < corpar.rows(); ++i)
    //     nll += pow(abs(corpar(i,j)), penalty) * lambda(2);
    // // nll += lambda(2) * pow(normCor, 1.0/(double)penalty);
    
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
	    norm += pow(fabs(efd(cc,i,j)), penalty); 
	  }
	  if(cc == 3)
	    norm += pow(fabs(efd(cc,0,j)), penalty); 
	  nll += lambda(0) * pow(norm, 1.0/(double)penalty);	  
	}
      }
    }
  }

  // Predict
  matrix<Type> logpred(NLEVELS(G),X_pred.cols());
  matrix<Type> logpredNoPrior(NLEVELS(G),X_pred.cols());
  for(int i = 0; i < logpred.cols(); ++i){
    vector<Type> lp(NLEVELS(G));
    Type lps = 0.0;
    for(int j = 0; j < NLEVELS(G); ++j){
      array<Type> tmpPred;
      tmpPred = X_pred.col(i);
      switch(modelType) {
      case Model::Independent:
        tmpPred -= (vector<Type>)(covar_pred.col(i).transpose() * muUse.col(j).matrix());
	if(commonCovar.rows() > 0){
	  tmpPred -= (vector<Type>)(commonCovar_pred.col(i).transpose() * commonMu);
	}
	break;
      case Model::Unstructured:
	tmpPred -= (vector<Type>)(covar_pred.col(i).transpose() * muUse.col(j).matrix());
	if(commonCovar.rows() > 0){
	  tmpPred -= (vector<Type>)(commonCovar_pred.col(i).transpose() * commonMu);
	}
	break;
      case Model::GMRF:
        tmpPred -= coords(j).array();
	break;
      default:
	Rf_error("Unknown model type");
	break;
      }
      logpredNoPrior(j,i) = -dist(j)(tmpPred);
      lp(j) = -dist(j)(tmpPred) + log(prior(j));
      if(j == 0){
	lps = lp(j);
      }else{
	lps = logspace_add(lps,lp(j));
      }
    }
    logpred.col(i) = lp - lps;
  }


  if(estimateUnbiasedTestProportions){
    matrix<Type> theta(thetaIn.rows() + 1,thetaIn.cols());
    for(int j = 0; j < theta.cols(); ++j){
      for(int i = 0; i < thetaIn.rows(); ++i)
	theta(i,j) = exp(thetaIn(i,j));
      theta(theta.rows()-1,j) = 1.0;
      theta.col(j) /= theta.col(j).sum();
    }
    REPORT(theta);
    ADREPORT(theta);

    if(trainContributionToTestProportions){
      for(int i = 0; i < G.size(); ++i)
	nll -= log(theta(G(i),trainBiasCorrectionGroup(i)));
    }

    matrix<Type> logpredBiasCorrected(NLEVELS(G),X_pred.cols());

    for(int i = 0; i < logpredNoPrior.cols(); ++i){
      vector<Type> th = theta.col(testBiasCorrectionGroup(i));
      vector<Type> lp(NLEVELS(G));
      Type lps = 0.0;
      for(int j = 0; j < NLEVELS(G); ++j){
	lp(j) = logpredNoPrior(j,i) + log(th(j));
	if(j == 0){
	  lps = lp(j);
	}else{
	  lps = logspace_add(lps,lp(j));
	}
      }
      if(addMisfitCategory)
	lps = logspace_add(lps,log(th(th.size()-1)));
      nll -= lps;
      logpredBiasCorrected.col(i) = lp - lps;
    }
    REPORT(logpredBiasCorrected);
  }

  
  REPORT(logpred);
  REPORT(logpredNoPrior);
  REPORT(coords);
  REPORT(meanCoords);
  REPORT(muUse);
  ADREPORT(muUse);
  
  return nll;
}
