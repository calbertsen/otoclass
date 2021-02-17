#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include "../inst/include/efd.hpp"

#define install Rf_install
#define findVar Rf_findVar
#define setAttrib Rf_setAttrib
#include <TMB.hpp>


namespace my_atomic {

  template<class T>
  T logspace_add2_raw (const T &logx, const T &logy) {
    // Was:
    //  fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
    if(logx == R_NegInf)
      return(logy);
    if(logy == R_NegInf)
      return(logx);
    return ( logx < logy ?
             logy + log1p (exp (logx - logy)) :
             logx + log1p (exp (logy - logx)) );
  }

  
TMB_BIND_ATOMIC(logspace_add2,
                11,
                logspace_add2_raw(x[0], x[1]) )

  
}

template<class Type>
Type logspace_add2(Type logx, Type logy) {
  if ( !CppAD::Variable(logx) && logx == Type(-INFINITY) )
    return logy;
  if ( !CppAD::Variable(logy) && logy == Type(-INFINITY) )
    return logx;
  CppAD::vector<Type> tx(3);
  tx[0] = logx;
  tx[1] = logy;
  tx[2] = 0; // order
  return my_atomic::logspace_add2(tx)[0];
}



template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x)) || asDouble(x) == NA_INTEGER;
}
bool isNA(int x){
  return x == NA_INTEGER;
}



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




//////////////////////////////////////////////////////////////////
/////////////// From fishfollower/stockassessment ///////////////
////////////////////////////////////////////////////////////////


template<class Type>
Type logspace_add_p (Type logx, Type logy, Type p) {
  return log((Type(1)-p)*exp(logy-logx)+p)+logx; // the order of x and y is taylored for this application 
}

template<class Type>
Type logdrobust(Type x, Type p, Type df){
  Type ld1=dnorm(x,Type(0.0),Type(1.0),true);
  if(p<Type(1.0e-16)){
    return ld1;
  }else{
    Type ld2=dt(x,df,true);
    Type logres=logspace_add_p(ld2,ld1,p);
    return logres;
  }
}
VECTORIZE3_ttt(logdrobust)

template <class Type>
class MVMIX_t{
  Type halfLogDetS;         
  vector<Type> p1;                  /*fraction t3*/
  vector<Type> df;
  matrix<Type> Sigma;       
  vector<Type> sd;
  matrix<Type> L_Sigma;
  matrix<Type> inv_L_Sigma;
public:
  MVMIX_t(){}
  MVMIX_t(matrix<Type> Sigma_, vector<Type> p1_, vector<Type> df_){
    setSigma(Sigma_);
    p1=p1_;
    df=df_;
  }
  matrix<Type> cov(){return Sigma;}
  void setSigma(matrix<Type> Sigma_){
    Sigma = Sigma_;
    sd = sqrt(vector<Type>(Sigma.diagonal()));
    Eigen::LLT<Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic> > llt(Sigma);
    L_Sigma = llt.matrixL();
    vector<Type> D=L_Sigma.diagonal();
    halfLogDetS = sum(log(D));
    inv_L_Sigma = atomic::matinv(L_Sigma);
  }
  void setSigma(matrix<Type> Sigma_, vector<Type> p1_){
    setSigma(Sigma_);
    p1=p1_;
  }
  /** \brief Evaluate the negative log density */
  Type operator()(vector<Type> x){
    vector<Type> z = inv_L_Sigma*x;
    return -sum(logdrobust(z,p1,df))+halfLogDetS;
  }
  // Type operator()(vector<Type> x, vector<Type> keep){
  //   matrix<Type> S = Sigma;
  //   vector<Type> not_keep = Type(1.0) - keep;
  //   for(int i = 0; i < S.rows(); i++){
  //     for(int j = 0; j < S.cols(); j++){
  // 	S(i,j) = S(i,j) * keep(i) * keep(j);
  //     }
  //     //S(i,i) += not_keep(i) * pow((Type(1)-p1)*sqrt(Type(0.5)/M_PI)+p1*(Type(1)/M_PI),2); //(t(1))
  //     S(i,i) += not_keep(i) * pow((Type(1)-p1)*sqrt(Type(0.5)/M_PI)+p1*(Type(2)/(M_PI*sqrt(df))),2);
  //   }
  //   return MVMIX_t<Type>(S,p1)(x * keep);
  // }

  vector<Type> simulate() {
    int siz = Sigma.rows();
    vector<Type> x(siz);
    for(int i=0; i<siz; ++i){
      Type u = runif(0.0,1.0);
      if(u<p1(i)){
        x(i) = rt(asDouble(df(i)));
      }else{
        x(i) = rnorm(0.0,1.0);
      }
    }
    x = L_Sigma*x;
    return x;
  }
};

template <class Type>
MVMIX_t<Type> MVMIX(matrix<Type> Sigma, vector<Type> p1){
  return MVMIX_t<Type>(Sigma,p1);
}


/////////////////////////////////////////////////////////////////////////////////////
/////////////////// End of code from fishfollower/stockassesment ///////////////////
///////////////////////////////////////////////////////////////////////////////////




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
  // SEPARABLE_t<GMRF_t<Type>, VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > > gmrf;

public:
  OTOLITH_t();
  OTOLITH_t(vector<Type> corpar,
	    vector<Type> sigma,
	    // SparseMatrix<Type> Q,
	    Model::Types model_);
  Type operator()(array<Type> x);

};

template<class Type>
OTOLITH_t<Type>::OTOLITH_t(){}

template<class Type>
OTOLITH_t<Type>::OTOLITH_t(vector<Type> corpar,
			   vector<Type> sigma,
			   // SparseMatrix<Type> Q,
			   Model::Types model_) :
  model(model_)
{
  switch(model) {
  case Model::Independent:
    independent = VECSCALE(N01<Type>(), sigma);
    break;
  case Model::Unstructured:
    unstructured = VECSCALE(UNSTRUCTURED_CORR(corpar), sigma);
    break;
  // case Model::GMRF:
  //   gmrf = SEPARABLE(GMRF(Q),VECSCALE(UNSTRUCTURED_CORR(corpar), sigma));
  //   break;
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
  // case Model::GMRF:
  //   r= gmrf(x);
  //   break;
  default:
    Rf_error("Unknown model type");
    break;
  }
  return r;
}


template <class Type>
Type objective_function<Type>::operator() () {

  DATA_INTEGER(model);
  Model::Types modelType = static_cast<Model::Types>(model);
    
  DATA_FACTOR(G);
  DATA_ARRAY(Y);
  DATA_SPARSE_MATRIX(X);
  DATA_SPARSE_MATRIX(XCom);
  DATA_FACTOR(proportionGroup);
  DATA_FACTOR(confusionGroup);
  DATA_FACTOR(dispersionGroup);
  
  DATA_INTEGER(penalty);

  DATA_ARRAY(Y_pred)
  DATA_SPARSE_MATRIX(X_pred);
  DATA_SPARSE_MATRIX(XCom_pred);
  DATA_FACTOR(proportionGroup_pred);
  DATA_FACTOR(dispersionGroup_pred);

  
  PARAMETER_ARRAY(mu);
  PARAMETER_MATRIX(commonMu);
  PARAMETER_MATRIX(logSigma);
  PARAMETER_MATRIX(corpar);

  PARAMETER_VECTOR(logLambda);

  PARAMETER_MATRIX(thetaIn);

  PARAMETER_ARRAY(MIn);

  PARAMETER_MATRIX(tmixpIn);
  PARAMETER_MATRIX(logDf);

  PARAMETER(logSdDispersion);
  PARAMETER_MATRIX(dispersion);

  ////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// Check //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  if(G.size() != Y.cols())
    Rf_error("Wrong input length. G and Y must match");

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// Prepare /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  
  // Transform lambda
  vector<Type> lambda = exp(-logLambda);

  // Transform df
  matrix<Type> df(logDf.rows(), logDf.cols());
    for(int i = 0; i < df.rows(); ++i)
      for(int j = 0; j < df.cols(); ++j)
	df(i,j) = exp(logDf(i,j));
  
  // Transform mu
  array<Type> muUse(mu.dim);
  for(int i = 0; i < muUse.dim[0]; ++i)
    for(int j = 0; j < muUse.dim[1]; ++j){
      muUse(i,j,0) = mu(i,j,0);
      for(int k = 1; k < muUse.dim[2]; ++k){
	muUse(i,j,k) = muUse(i,j,0) + mu(i,j,k);
      }
    }
    
  // Transform sigma
  matrix<Type> sigma(logSigma.rows(), logSigma.cols());
  for(int i = 0; i < sigma.rows(); ++i)
    for(int j = 0; j < sigma.cols(); ++j)
      sigma(i,j) = exp(logSigma(i,j));

  
  // Construct covariance matrices
  vector<matrix<Type> > SigmaList(logSigma.cols());
  for(int k = 0; k < sigma.cols(); ++k){
    matrix<Type> stmp(sigma.rows(),sigma.rows());
    stmp.setIdentity();
    if(stmp.cols() > 1){
      stmp = UNSTRUCTURED_CORR((vector<Type>)corpar.col(k)).cov();
    }
    for(int i = 0; i < sigma.rows(); ++i)
      for(int j = 0; j < sigma.rows(); ++j)
	stmp(i,j) *= sigma(i,k) * sigma(j,k);
    SigmaList(k) = stmp;
  }


  
  // Transform theta
  matrix<Type> theta(thetaIn.rows() + 1,thetaIn.cols());
  for(int k = 0; k < dispersion.rows(); ++k){
    for(int j = 0; j < theta.cols(); ++j){
      for(int i = 0; i < thetaIn.rows(); ++i)
	theta(i,j) = exp(thetaIn(i,j));
      theta(theta.rows()-1,j) = 1.0;
      theta.col(j) /= theta.col(j).sum();
    }
  }

  // Transform confusion
  vector<matrix<Type> > Mvec(MIn.dim[2]);
  for(int k = 0; k < MIn.dim[2]; ++k){
    matrix<Type> Mtmp(MIn.dim[0],MIn.dim[1]);
    Mtmp.setZero();
    for(int j = 0; j < Mtmp.cols(); ++j){
      for(int i = 0; i < Mtmp.rows(); ++i){
	if(i != j){
	  Mtmp(i,j) = exp(MIn(i,j,k));
	  Mtmp(j,j) -= Mtmp(i,j);
	}
      }
    }
    Mvec(k) = expm(Mtmp);
  }

  // Student's t mixture
  matrix<Type> tmixp(tmixpIn.rows(),tmixpIn.cols());
  for(int i = 0; i < tmixpIn.rows(); ++i)
    for(int j = 0; j < tmixpIn.cols(); ++j)
      tmixp(i,j) = 1.0 / (1.0 + exp(-tmixpIn(i,j)));
  
  // Prepare observational distributions
  vector<MVMIX_t<Type> > dist(NLEVELS(G));
  for(int i = 0; i < NLEVELS(G); ++i)
    dist(i) = MVMIX_t<Type>(SigmaList(i), (vector<Type>)tmixp.col(i), (vector<Type>)df.col(i));
    //    dist(i) = OTOLITH_t<Type>((vector<Type>)corpar.col(i), (vector<Type>)sigma.col(i), modelType);

  Type sdDispersion = exp(logSdDispersion);
  
  // Likelihood value
  Type nll = 0.0;

  // Dispersion random effect
    for(int i = 0; i < dispersion.rows(); ++i){
      for(int j = 0; j < dispersion.cols(); ++j){
	if(CppAD::Variable(dispersion(i,j)))
	  nll -= dnorm(dispersion(i,j), Type(0.0), sdDispersion, true);      
      }
    }
  
  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Estimate /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  array<Type> posterior_mean(Y.rows(), Y.cols(), NLEVELS(G));
  posterior_mean.setZero();
  matrix<Type> prior_logprobability(NLEVELS(G), Y.cols());
  prior_logprobability.setZero();
  matrix<Type> posterior_logprobability(NLEVELS(G), Y.cols());
  posterior_logprobability.setZero();
  
  for(int i = 0; i < Y.cols(); ++i){ // Loop over individuals
    //vector<Type> th = thetaDisp.col(dispersionGroup(i)).col(proportionGroup(i));
    vector<Type> th(NLEVELS(G));
    th.setZero();
    if(isNA(dispersionGroup(i))){
      th = theta.col(proportionGroup(i));
    }else{
      for(int qq = 0; qq < thetaIn.rows(); ++qq)
	th(qq) = exp(thetaIn(qq,proportionGroup(i)) + dispersion(qq,dispersionGroup(i)));
      th(th.size() - 1) = 1.0;
      th /= th.sum();
    }
    //prior_logprobability.col(i) = log(th);
    vector<Type> Muse = Mvec(confusionGroup(i)).row(G(i));
    Type NormalizationConstant = (Muse*th).sum();
    vector<Type> lp(NLEVELS(G));
    Type lps = R_NegInf;
    Type postProbNorm = R_NegInf;
    // Contribution from G
    nll -= log(NormalizationConstant);

    // Contrbution from Y
    for(int j = 0; j < NLEVELS(G); ++j){ // Loop over groups
      array<Type> tmp;
      tmp = Y.col(i);
      vector<Type> tmpMean = (vector<Type>)(X.col(i).transpose() * muUse.col(j).matrix());
      if(XCom.rows() > 0){
	tmpMean += (vector<Type>)(XCom.col(i).transpose() * commonMu);
      }
      tmp -= tmpMean;
      posterior_mean.col(j).col(i) = tmpMean;
      Type ld = -dist(j)(tmp); // dist returns negative log-likelihood!
      lp(j) = ld + log(th(j)) + log(Muse(j)) - log(NormalizationConstant);
      posterior_logprobability(j,i) = ld + log(th(j));
      lps = logspace_add2(lps,lp(j));
      postProbNorm = logspace_add2(postProbNorm, posterior_logprobability(j,i));      
    }
    for(int j = 0; j < NLEVELS(G); ++j)
      posterior_logprobability(j,i) -= postProbNorm;
    nll -= lps;
  }

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Regularization ///////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  if(penalty != 0){
    for(int i = 0; i < mu.dim[0]; ++i)
      for(int j = 0; j < mu.dim[1]; ++j){
	for(int k = 0; k < mu.dim[2]; ++k){
	  switch(penalty){
	  case -1:
	    nll -= dt(mu(i,j,k) / lambda(0),Type(3.0),true) - log(lambda(0));
	    break;
	  case 1:
	    nll -= -fabs(mu(i,j,k)) / sigma(j,k) / lambda(0) - log(2.0 * lambda(0));
	    break;
	  case 2:
	    nll -= dnorm(mu(i,j,k),Type(0.0),lambda(0) * sigma(j,k), true);
	    break;
	  default:
	    nll -= my_atomic::Lp(mu(i,j,k),(Type)penalty,lambda(0) * sigma(j,k));
	  }
	}
      }
    vector<Type> sigmaMeans = sigma.rowwise().mean();
    for(int i = 0; i < commonMu.rows(); ++i)
      for(int j = 0; j < commonMu.cols(); ++j){
	switch(penalty){
	case -1:
	  nll -= dt(commonMu(i,j) / lambda(1),Type(3.0),true) - log(lambda(1));
	  break;
	case 1:
	  nll -= -fabs(commonMu(i,j)) / sigmaMeans(j) / lambda(1) - log(2.0 * lambda(1));
	  break;
	case 2:
	  nll -= dnorm(commonMu(i,j),Type(0.0),lambda(1) * sigmaMeans(j), true);
	  break;
	default:
	  nll -= my_atomic::Lp(commonMu(i,j),(Type)penalty,lambda(1) * sigmaMeans(j));
	}
      }
  }

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Predict //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  array<Type> pred_posterior_mean(Y_pred.rows(), Y_pred.cols(), NLEVELS(G));
  matrix<Type> pred_prior_logprobability(NLEVELS(G), Y_pred.cols());
  matrix<Type> pred_posterior_logprobability(NLEVELS(G), Y_pred.cols());
  
  for(int i = 0; i < Y_pred.cols(); ++i){ // Loop over individuals
    vector<Type> th(thetaIn.rows()+1);
    th.setZero();;
    if(isNA(dispersionGroup_pred(i))){
      th = theta.col(proportionGroup_pred(i));
    }else{
      for(int qq = 0; qq < thetaIn.rows(); ++qq)
	th(qq) = exp(thetaIn(qq,proportionGroup_pred(i)) + dispersion(qq,dispersionGroup_pred(i)));
      th(th.size() - 1) = 1.0;
      th /= th.sum();
    }
    Type postProbNorm = R_NegInf;
    pred_prior_logprobability .col(i) = th;
    for(int j = 0; j < NLEVELS(G); ++j){ // Loop over groups
      array<Type> tmp;
      tmp = Y_pred.col(i);
      vector<Type> tmpMean = (vector<Type>)(X_pred.col(i).transpose() * muUse.col(j).matrix());
      if(XCom_pred.rows() > 0){
	tmpMean += (vector<Type>)(XCom_pred.col(i).transpose() * commonMu);
      }
      tmp -= tmpMean;
      pred_posterior_mean.col(j).col(i) = tmpMean;
      Type ld =  dist(j)(tmp);
      pred_posterior_logprobability(j,i) = ld + log(th(j));
      postProbNorm = logspace_add2(postProbNorm, pred_posterior_logprobability(j,i));
    }
    for(int j = 0; j < NLEVELS(G); ++j)
      pred_posterior_logprobability(j,i) -= postProbNorm;
  }
  
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// Output //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  REPORT(muUse);
  ADREPORT(muUse);
  ADREPORT(commonMu);

  REPORT(posterior_mean);
  REPORT(posterior_logprobability);
  REPORT(prior_logprobability);

  ADREPORT(posterior_mean);
  ADREPORT(posterior_logprobability);

  REPORT(pred_posterior_mean);
  REPORT(pred_posterior_logprobability);
  REPORT(pred_prior_logprobability);

  ADREPORT(pred_posterior_mean);
  ADREPORT(pred_posterior_logprobability);


  REPORT(theta);
  ADREPORT(theta);
  
  REPORT(Mvec);

  REPORT(SigmaList);

  REPORT(tmixp);

  REPORT(df);
  
  return nll;
  }
