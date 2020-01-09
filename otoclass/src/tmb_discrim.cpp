#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include "../inst/include/efd.hpp"

#define install Rf_install
#define findVar Rf_findVar
#define setAttrib Rf_setAttrib
#define warning Rf_warning
#define R_NO_REMAP
#include <TMB.hpp>
#include "../inst/include/Rmath/toTMB.hpp"

#include <fenv.h>

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
}

template<class Type>
Type Lp(Type x, Type p, Type lambda) {
  vector<Type> args(4); // Last index reserved for derivative order
  args << x, p, lambda, 0;
  return my_atomic::func(CppAD::vector<Type>(args))[0];
}


//////////////////////////////////////////////////////////////////
/////////////// From fishfollower/stockassessment ///////////////
////////////////////////////////////////////////////////////////



// template<class Type>
// Type logspace_add_p (Type logx, Type logy, Type p) {
//   return log((Type(1.0)-p)*exp(logy-logx)+p)+logx; // the order of x and y is taylored for this application
// }

// template<class Type>
// Type logdrobust(Type x, Type p, Type df){
//   Type ld1=dnorm(x,Type(0.0),Type(1.0),true);
//   Type ld2=dt(x,df,true);
//   Type logres=logspace_add_p(ld2,ld1,p);
//   return logres;
//   // return CppAD::CondExpLt(p,Type(1e-8),ld1,
//   // 			  CppAD::CondExpGt(p,Type(0.99999999),ld2,logres));
// }
// VECTORIZE3_ttt(logdrobust)


// template<class Type>
// Type logprobust(Type x, Type p, Type df){
//   Type ld1=pnorm5(x,Type(0.0),Type(1.0),1,1);
//   Type ld2=pt(x,df,1,1);
//   Type logres=logspace_add_p(ld2,ld1,p);
//   return logres;
//   // return CppAD::CondExpLt(p,Type(1e-8),ld1,
//   // 			  CppAD::CondExpGt(p,Type(0.99999999),ld2,logres));
// }
// VECTORIZE3_ttt(logprobust)



namespace my_atomic {

  
  template<class Float>
  Float logspace_add_p_raw (Float logx, Float logy, Float logit_p) {
    Float NINF = ML_NEGINF;
    Float PINF = ML_POSINF;
    Float zero = 0;
    // if(logit_p < -7)
    //   logit_p = -7;
    // if(logit_p < 7)
    //   logit_p = 7;
    Float log_p   = - atomic::robust_utils::logspace_add( zero , Float(-logit_p) );
    Float log_1mp = - atomic::robust_utils::logspace_add( zero , logit_p );
    if(log_p == NINF) // p == 0
      return logy;
    if(log_1mp == PINF) // 1-p == 0 <=> p == 1
      return logx;
    
    Float v1 = logx + log_p; 
    Float v2 = logy + log_1mp;
    
    if(!R_FINITE(v1) && !R_FINITE(v2)){
      return NINF;
    }else if(R_FINITE(v1) && !R_FINITE(v2)){
      return v1;
    }else if(!R_FINITE(v1) && R_FINITE(v2)){
      return v2;
    }else{      
      Float finalV = atomic::robust_utils::logspace_add(v1, v2);
      return finalV;   
    }
  }
  TMB_BIND_ATOMIC(logspace_add_p1, 111, logspace_add_p_raw(x[0], x[1], x[2]))


  template<class Float>
  Float dkwnorm_raw(Float x, Float a, Float b){
    Float one = 1;
    Float zero = 0;
    Float eps = std::numeric_limits<double>::epsilon();
    Float logg = dnorm(x, zero, one, true);
    Float logG = pnorm5_raw(x, zero, one, one, one);
    if(!R_FINITE(logg)) // if logg is -inf, dkwmix is -inf
      return logg;
    if(!R_FINITE(logG)) // if logG is -inf, dkwmix is -inf
      return logG;
    Float v0 = a*logG;
    Float v1 = atomic::robust_utils::logspace_sub(zero, Float(v0 - eps));
    if(!R_FINITE(v1)) // if (a * logG == 0), dkwmix is -inf
      return v1;
    Float v = log(a) + log(b) + logg + (a-one) * logG + (b-one) * v1;
    return v;
  }

  template<class Float>
  Float pkwnorm_raw(Float x, Float a, Float b){
    Float NINF = ML_NEGINF;
    Float eps = std::numeric_limits<double>::epsilon();
    Float one = 1;
    Float zero = 0;
    // Float logg = dnorm(x,zero,one,true);
    Float logG = pnorm5_raw(x,zero,one,one,one);
    if(logG == NINF)
      return logG;
    Float v1 = a * logG;
    Float v2 = b * atomic::robust_utils::logspace_sub(zero, Float(v1 - eps));
    if(!R_FINITE(v2)) // log(1-exp(a*log(G(x)))) <=> G(x)^a == 1 => F(x) == 1 <=> log(F(x) == 0; logspace_sub can return NaN for logspace_sub(0, 0) if the second 0 is slightly larger
      return zero;
    Float v3 = atomic::robust_utils::logspace_sub(zero,Float(v2 - eps));
    if(!R_FINITE(v3)) // log(1-exp(a*log(G(x)))) <=> G(x)^a == 1 => F(x) == 1 <=> log(F(x) == 0; logspace_sub can return NaN for logspace_sub(0, 0) if the second 0 is slightly larger
      return zero;
    return v3;
  }

  
  template<class Float>
  Float logdrobust_raw(Float x, Float logit_p, Float df, Float a, Float b){
    Float ld1=dkwnorm_raw(x,a,b);
    Float ld2=dt(x,df,true);
    Float logres=logspace_add_p_raw(ld2,ld1,logit_p);
    return logres;
  }

  TMB_BIND_ATOMIC(logdrobust1, 11111, logdrobust_raw(x[0], x[1], x[2], x[3], x[4]))

  template<class Float>
  Float logprobust_raw(Float x, Float logit_p, Float df, Float a, Float b){
    Float ld1=pkwnorm_raw(x,a,b);
    Float ld2=pt_raw(x,df,Float(1),Float(1));
    Float logres=logspace_add_p_raw(ld2,ld1,logit_p);
    return logres;
  }

  TMB_BIND_ATOMIC(logprobust1, 11111, logprobust_raw(x[0], x[1], x[2], x[3], x[4]))
  
  // template<class Float>
  // Float dkwmix_raw(Float x, Float logit_p, Float df, Float a, Float b){
  //   Float logg = logdrobust_raw(x,logit_p,df);
  //   Float logG = logprobust_raw(x,logit_p,df);
  //   if(!R_FINITE(logg)) // if logg is -inf, dkwmix is -inf
  //     return logg;
  //   if(!R_FINITE(logG)) // if logG is -inf, dkwmix is -inf
  //     return logG;
  //   if(fabs(logG) < 1e-8 || a < 1e-8) // if a * logG == 0, dkwmix is -inf
  //     return ML_NEGINF;
  //   Float v = log(a) + log(b) + logg + (a-Float(1.0)) * logG + (b-Float(1.0)) * atomic::robust_utils::logspace_sub(Float(0.0), Float(a*logG));
  //   return v;
  // }
  // TMB_BIND_ATOMIC(dkwmix1, 11111, dkwmix_raw(x[0], x[1], x[2], x[3], x[4]))


  // template<class Float>
  // Float pkwmix_raw(Float x, Float logit_p, Float df, Float a, Float b){
  //   Float logG = logprobust_raw(x,logit_p,df);
  //   // Float v = log(1.0 - pow(1.0 - exp(a * logG), b));
  //   // return v;
  //   Float v1 = a * logG;
  //   if(!R_Finite(log(1-exp(v1))))
  //     return Float(0.0);
  //   Float v2 = b * atomic::robust_utils::logspace_sub(Type(0.0), v1);
  //   if(!R_Finite(log(1-exp(v2))))
  //     return Float(0.0);
  //   return atomic::robust_utils::logspace_sub(Float(0.0),v2);
  //   // Type v =  atomic::robust_utils::logspace_sub(Float(0.0), b * atomic::robust_utils::logspace_sub(Type(0.0), a * logG));
  //   // if(ISNAN(v))
  //   //   v = 0.0;
  //   // return v;
  // }
  // TMB_BIND_ATOMIC(pkwmix1, 11111, pkwmix_raw(x[0], x[1], x[2], x[3], x[4]))

 
  
  template<class Float>
  struct KWMOMENT_t {
    typedef Float Scalar;
    Float a;
    Float b;
    int m;
    Float operator() (Float x) {
      Float v = exp(dkwnorm_raw(x,a,b));
      Float v2 =  pow(x, m) * v;
      if(ISNAN(v))
	v = Float(0);
      return v2;
    }
    Float get() {
      using gauss_kronrod::integrate;
      Float ans = integrate(*this, -INFINITY, INFINITY);
      return ans;
    }
  };

  template<class Float>
  Float kwmoment_raw(Float a, Float b, Float m) {
    KWMOMENT_t<Float> f = {a,b, (int) trunc(m)};
    Float v = f.get();
    return v;
  }
  TMB_BIND_ATOMIC(kwmoment1, 110, kwmoment_raw(x[0], x[1], x[2]))


   
}


template<class Type>
Type logspace_add_p(Type logx, Type logy, Type logit_p){
  vector<Type> args(4); // Last index reserved for derivative order
  args << logx, logy, logit_p, 0;
  return my_atomic::logspace_add_p1(CppAD::vector<Type>(args))[0];
}


template<class Type>
Type logdkwmix(Type x, Type logit_p, Type df, Type a, Type b){
  vector<Type> args(6); // Last index reserved for derivative order
  args << x, logit_p, df, a, b, 0;
  return my_atomic::logdrobust1(CppAD::vector<Type>(args))[0];
}

GVECTORIZE(logdkwmix,V,V,V,V,V,N)

template<class Type>
Type logpkwmix(Type x, Type logit_p, Type df, Type a, Type b){
  vector<Type> args(6); // Last index reserved for derivative order
  args << x, logit_p, df, a, b, 0;
  return my_atomic::logprobust1(CppAD::vector<Type>(args))[0];
}

GVECTORIZE(logpkwmix,V,V,V,V,V,N)


template<class Type>
Type kwmoment(Type a, Type b, int m) {
  vector<Type> args(4); // Last index reserved for derivative order
  args << a, b, m, 0;
  return my_atomic::kwmoment1(CppAD::vector<Type>(args))[0];
}

GVECTORIZE(kwmoment,V,V,I,N, N, N)



template <class Type>
class MVMIX_t{
  Type halfLogDetS;         
  vector<Type> logit_p;                  /*fraction t3*/
  vector<Type> df;
  vector<Type> a;
  vector<Type> b;
  matrix<Type> Sigma;       
  vector<Type> sd;
  matrix<Type> L_Sigma;
  matrix<Type> inv_L_Sigma;
  vector<Type> KWm1;
  vector<Type> KWm2;
  vector<Type> KWsd;

  void getKWMoments(){
    KWm1 = kwmoment(a, b, 1);
    KWm2 = kwmoment(a, b, 2);
    KWsd = sqrt(KWm2 - KWm1 * KWm1);    
  }

public:
  MVMIX_t(){}
  MVMIX_t(matrix<Type> Sigma_, vector<Type> logit_p_, vector<Type> df_, vector<Type> a_, vector<Type> b_) :
    halfLogDetS(0), logit_p(logit_p_), df(df_), a(a_), b(b_), Sigma(Sigma_), sd(Sigma_.cols()), L_Sigma(Sigma_), inv_L_Sigma(Sigma_), KWm1(logit_p_.size()), KWm2(logit_p_.size()), KWsd(logit_p_.size()) {
    setSigma(Sigma_);
    getKWMoments();
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
  /** \brief Evaluate the negative log density */
  Type operator()(vector<Type> x){
    vector<Type> z = inv_L_Sigma*x;
    //return -sum(logdrobust(z,p1,df))+halfLogDetS;
    vector<Type> z2 = z * KWsd + KWm1;
    return -sum(logdkwmix(z2,logit_p,df,a,b)) + halfLogDetS - sum(log(KWsd));
    
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

  // vector<Type> simulate() {
  //   int siz = Sigma.rows();
  //   vector<Type> x(siz);
  //   for(int i=0; i<siz; ++i){
  //     Type u = runif(0.0,1.0);
  //     if(u<p1(i)){
  //       x(i) = rt(asDouble(df(i)));
  //     }else{
  //       x(i) = rnorm(0.0,1.0);
  //     }
  //   }
  //   x = L_Sigma*x;
  //   return x;
  // }
};


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


// template<class Type>
// class OTOLITH_t {
  
// private:
//   Model::Types model;
//   VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > unstructured;
//   VECSCALE_t<N01<Type> > independent;
//   // SEPARABLE_t<GMRF_t<Type>, VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > > gmrf;

// public:
//   OTOLITH_t();
//   OTOLITH_t(vector<Type> corpar,
// 	    vector<Type> sigma,
// 	    // SparseMatrix<Type> Q,
// 	    Model::Types model_);
//   Type operator()(array<Type> x);

// };

// template<class Type>
// OTOLITH_t<Type>::OTOLITH_t(){}

// template<class Type>
// OTOLITH_t<Type>::OTOLITH_t(vector<Type> corpar,
// 			   vector<Type> sigma,
// 			   // SparseMatrix<Type> Q,
// 			   Model::Types model_) :
//   model(model_)
// {
//   switch(model) {
//   case Model::Independent:
//     independent = VECSCALE(N01<Type>(), sigma);
//     break;
//   case Model::Unstructured:
//     unstructured = VECSCALE(UNSTRUCTURED_CORR(corpar), sigma);
//     break;
//   // case Model::GMRF:
//   //   gmrf = SEPARABLE(GMRF(Q),VECSCALE(UNSTRUCTURED_CORR(corpar), sigma));
//   //   break;
//   default:
//     Rf_error("Unknown model type");
//     break;
//   }
// }

// template<class Type>
// Type OTOLITH_t<Type>::operator() (array<Type> x){
//   Type r = 0.0;
//   switch(model) {
//   case Model::Independent:
//     r = independent((vector<Type>)x);
//     break;
//   case Model::Unstructured:
//     r= unstructured((vector<Type>)x);
//     break;
//   // case Model::GMRF:
//   //   r= gmrf(x);
//   //   break;
//   default:
//     Rf_error("Unknown model type");
//     break;
//   }
//   return r;
// }


template <class Type>
Type objective_function<Type>::operator() () {

  DATA_INTEGER(model);
  Model::Types modelType = static_cast<Model::Types>(model);
    
  DATA_FACTOR(G);
  DATA_MATRIX(Y);
  DATA_SPARSE_MATRIX(X);
  DATA_SPARSE_MATRIX(XCom);
  DATA_FACTOR(proportionGroup);
  DATA_FACTOR(confusionGroup);

  DATA_INTEGER(penalty);

  DATA_MATRIX(Y_pred);
  DATA_SPARSE_MATRIX(X_pred);
  DATA_SPARSE_MATRIX(XCom_pred);
  DATA_FACTOR(proportionGroup_pred);

  DATA_INTEGER(kw_idparam);
  
  PARAMETER_ARRAY(mu);
  PARAMETER_MATRIX(commonMu);
  PARAMETER_MATRIX(logSigma);
  PARAMETER_MATRIX(corpar);

  PARAMETER_VECTOR(logLambda);

  PARAMETER_MATRIX(thetaIn);

  PARAMETER_ARRAY(MIn);

  PARAMETER_MATRIX(tmixpIn);
  PARAMETER_MATRIX(logDf);

  PARAMETER_MATRIX(kw1);
  PARAMETER_MATRIX(kw2);

  

  ////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// Check //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  if(G.size() != Y.cols())
    Rf_error("Wrong input length. Groups and Observations must match");
  if(Y.cols() != X.cols())
    Rf_error("Wrong input length. Observation and Design Matrix must match");
  if(mu.dim[0] != X.rows())
    Rf_error("Wrong input length. Design Matrix must match the number of mean parameters");
  if(mu.dim[1] != Y.rows())
    Rf_error("Wrong input length. Observations must match the number of mean parameters");
  if(mu.dim[2] != NLEVELS(G))
    Rf_error("Wrong input length. Groups and mean parameters must match");

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

  // Transform a
  matrix<Type> aParam(kw1.rows(), kw1.cols());
  for(int i = 0; i < aParam.rows(); ++i)
    for(int j = 0; j < aParam.cols(); ++j){
      aParam(i,j) = exp(kw1(i,j) + Type(kw_idparam) * (kw2(i,j)));
      //aParam(i,j) = exp(logA(i,j));
    }

  // Transform b
  matrix<Type> bParam(kw2.rows(), kw2.cols());
  for(int i = 0; i < bParam.rows(); ++i)
    for(int j = 0; j < bParam.cols(); ++j){
      bParam(i,j) = exp(kw2(i,j) + Type(kw_idparam) * (kw1(i,j) - 2.0 * kw2(i,j)));
      // bParam(i,j) = exp(logB(i,j));
    }
  
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
  for(int j = 0; j < theta.cols(); ++j){
    for(int i = 0; i < thetaIn.rows(); ++i)
      theta(i,j) = exp(thetaIn(i,j));
    theta(theta.rows()-1,j) = 1.0;
    theta.col(j) /= theta.col(j).sum();
  }

  // Transform confusion
  vector<matrix<Type> > Mvec(MIn.dim[2]);
  for(int k = 0; k < MIn.dim[2]; ++k){
    matrix<Type> Mtmp(MIn.dim[0],MIn.dim[1]);
    Mtmp.setZero();
    for(int j = 0; j < Mtmp.cols(); ++j){
      for(int i = 0; i < Mtmp.rows(); ++i){
	// if(i != j){
	//   Mtmp(i,j) = exp(MIn(i,j,k));
	//   Mtmp(j,j) -= Mtmp(i,j);
	// }
	Mtmp(i,j) = exp(MIn(i,j,k));
      }
      Mtmp.col(j) /= Mtmp.col(j).sum();
    }
    Mvec(k) = Mtmp; //expm(Mtmp);
  }

  // Student's t mixture
  matrix<Type> tmixp(tmixpIn.rows(),tmixpIn.cols());
  for(int i = 0; i < tmixpIn.rows(); ++i)
    for(int j = 0; j < tmixpIn.cols(); ++j)
      tmixp(i,j) = 1.0 / (1.0 + exp(-tmixpIn(i,j)));
  
  // Prepare observational distributions
  vector<MVMIX_t<Type> > dist(NLEVELS(G));
  for(int i = 0; i < NLEVELS(G); ++i)
    dist(i) = MVMIX_t<Type>(SigmaList(i),
			    (vector<Type>)tmixpIn.col(i),
			    (vector<Type>)df.col(i),
			    (vector<Type>)aParam.col(i),
			    (vector<Type>)bParam.col(i));
  //    dist(i) = OTOLITH_t<Type>((vector<Type>)corpar.col(i), (vector<Type>)sigma.col(i), modelType);

  
  // Likelihood value
  //Type nll = 0.0;
  parallel_accumulator<Type> nll(this);

  
  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Estimate /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  array<Type> posterior_mean(Y.rows(), Y.cols(), NLEVELS(G));
  matrix<Type> posterior_logprobability(NLEVELS(G), Y.cols());
  
  for(int i = 0; i < Y.cols(); ++i){ // Loop over individuals
    Type nll_est = 0.0;
    vector<Type> th = theta.col(proportionGroup(i));
    vector<Type> Muse = Mvec(confusionGroup(i)).row(G(i));
    Type NormalizationConstant = (Muse*th).sum();
    vector<Type> lp(NLEVELS(G));
    lp.setZero();
    Type lps = 0.0;
    Type postProbNorm = 0.0;
    // Contribution from G
    nll_est -= log(NormalizationConstant);
    // Contrbution from Y
    for(int j = 0; j < NLEVELS(G); ++j){ // Loop over groups
      vector<Type> tmp((vector<Type>)Y.col(i));
      vector<Type> tmpMean((vector<Type>)(X.col(i).transpose() * muUse.col(j).matrix()));
      if(XCom.rows() > 0){
	tmpMean += (vector<Type>)(XCom.col(i).transpose() * commonMu);
      }
      tmp -= tmpMean;
      posterior_mean.col(j).col(i) = tmpMean;
      Type ld = -dist(j)((vector<Type>)tmp); // dist returns negative log-likelihood!
      lp(j) = ld + log(squeeze(th(j))) + log(squeeze(Muse(j))) - log(squeeze(NormalizationConstant));
      posterior_logprobability(j,i) = ld + log(squeeze(th(j)));
      if(j == 0){
	lps = lp(j);
	postProbNorm = posterior_logprobability(j,i);
      }else{
	lps = logspace_add(lps,lp(j));
	// Avoid NaN for the -Inf + -Inf case
	//lps = CppAD::CondExpLt(lp(j),Type(-500.0),lps,logspace_add(lps,lp(j)));
	postProbNorm = logspace_add(postProbNorm, posterior_logprobability(j,i));
      }
    }
    for(int j = 0; j < NLEVELS(G); ++j)
      posterior_logprobability(j,i) -= postProbNorm;
    nll_est -= lps;
    nll += nll_est;
  }

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Regularization ///////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  if(penalty != 0){
    Type nll_pen = 0.0;
    for(int i = 0; i < mu.dim[0]; ++i)
      for(int j = 0; j < mu.dim[1]; ++j){
	for(int k = 0; k < mu.dim[2]; ++k){
	  switch(penalty){
	  case -1:
	    nll_pen -= dt(mu(i,j,k) / lambda(0),Type(3.0),true) - log(lambda(0));
	    break;
	  case 1:
	    nll_pen -= -fabs(mu(i,j,k)) / sigma(j,k) / lambda(0) - log(2.0 * lambda(0));
	    break;
	  case 2:
	    nll_pen -= dnorm(mu(i,j,k),Type(0.0),lambda(0) * sigma(j,k), true);
	    break;
	  default:
	    nll_pen -= Lp(mu(i,j,k),(Type)penalty,lambda(0) * sigma(j,k));
	  }
	}
      }
    vector<Type> sigmaMeans = sigma.rowwise().mean();
    for(int i = 0; i < commonMu.rows(); ++i)
      for(int j = 0; j < commonMu.cols(); ++j){
	switch(penalty){
	case -1:
	  nll_pen -= dt(commonMu(i,j) / lambda(1),Type(3.0),true) - log(lambda(1));
	  break;
	case 1:
	  nll_pen -= -fabs(commonMu(i,j)) / sigmaMeans(j) / lambda(1) - log(2.0 * lambda(1));
	  break;
	case 2:
	  nll_pen -= dnorm(commonMu(i,j),Type(0.0),lambda(1) * sigmaMeans(j), true);
	  break;
	default:
	  nll_pen -= Lp(commonMu(i,j),(Type)penalty,lambda(1) * sigmaMeans(j));
	}
      }
      nll += nll_pen;
  }

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Predict //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  array<Type> pred_posterior_mean(Y_pred.rows(), Y_pred.cols(), NLEVELS(G));
  matrix<Type> pred_posterior_logprobability(NLEVELS(G), Y_pred.cols());
  
  if( isDouble<Type>::value && TMB_OBJECTIVE_PTR -> current_parallel_region<0 ){
  
  for(int i = 0; i < Y_pred.cols(); ++i){ // Loop over individuals
    vector<Type> th = theta.col(proportionGroup_pred(i));
    Type postProbNorm = 0.0;

    for(int j = 0; j < NLEVELS(G); ++j){ // Loop over groups
      vector<Type> tmp((vector<Type>)Y_pred.col(i));
      vector<Type> tmpMean((vector<Type>)(X_pred.col(i).transpose() * muUse.col(j).matrix()));
      if(XCom_pred.rows() > 0){
	tmpMean += (vector<Type>)(XCom_pred.col(i).transpose() * commonMu);
      }
      tmp -= tmpMean;
      pred_posterior_mean.col(j).col(i) = tmpMean;
      Type ld =  dist(j)(tmp);
      pred_posterior_logprobability(j,i) = ld + log(th(j));
      if(j == 0){
	postProbNorm = pred_posterior_logprobability(j,i);
      }else{
	postProbNorm = logspace_add
	  (postProbNorm, pred_posterior_logprobability(j,i));
      }
    }
    for(int j = 0; j < NLEVELS(G); ++j)
      pred_posterior_logprobability(j,i) -= postProbNorm;
  }
  
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// Output //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  REPORT(muUse);
  ADREPORT(muUse);
  ADREPORT(commonMu);

  REPORT(posterior_mean);
  REPORT(posterior_logprobability);

  // ADREPORT(posterior_mean);
  // ADREPORT(posterior_logprobability);

  REPORT(pred_posterior_mean);
  REPORT(pred_posterior_logprobability);

  ADREPORT(pred_posterior_mean);
  // ADREPORT(pred_posterior_logprobability);


  REPORT(theta);
  ADREPORT(theta);

  REPORT(Mvec);

  REPORT(SigmaList);

  REPORT(tmixp);

  REPORT(df);

  REPORT(aParam);

  REPORT(bParam);
  
  return nll;
}





template<class Type>
Type logdkwmixE0(Type x, Type logit_p, Type df, Type a, Type b){
  Type KWm1 = kwmoment(a, b, 1);
  Type KWm2 = kwmoment(a, b, 2);
  Type KWsd = sqrt(KWm2 - KWm1 * KWm1);      
  Type y = x * KWsd + KWm1; // x should have mean zero and variance 1
  return logdkwmix(y,logit_p,df,a,b) + log(KWsd);
}

GVECTORIZE(logdkwmixE0,V,V,V,V,V,N)

template<class Type>
Type logpkwmixE0(Type x, Type logit_p, Type df, Type a, Type b){
  Type KWm1 = kwmoment(a, b, 1);
  Type KWm2 = kwmoment(a, b, 2);
  Type KWsd = sqrt(KWm2 - KWm1 * KWm1);      
  Type y = x * KWsd + KWm1; // x should have mean zero and variance 1
  return logpkwmix(y,logit_p,df,a,b);
}

GVECTORIZE(logpkwmixE0,V,V,V,V,V,N)




extern "C" {
  SEXP logspace_add_p (SEXP logx, SEXP logy, SEXP logit_p){
    SEXP res;
    PROTECT(res = Rf_allocVector(REALSXP,1));
    REAL(res)[0] = logspace_add_p<double>(Rf_asReal(logx), Rf_asReal(logy), Rf_asReal(logit_p));
    UNPROTECT(1);
    return res;
  }
  
  SEXP logdkwmix(SEXP x, SEXP logit_p, SEXP df, SEXP a, SEXP b){
    SEXP res;
    PROTECT(res = Rf_allocVector(REALSXP,1));
    REAL(res)[0] = logdkwmix<double>(Rf_asReal(x), Rf_asReal(logit_p), Rf_asReal(df), Rf_asReal(a), Rf_asReal(b));
    UNPROTECT(1);
    return res;
  }
  SEXP logpkwmix(SEXP x, SEXP logit_p, SEXP df, SEXP a, SEXP b){
    SEXP res;
    PROTECT(res = Rf_allocVector(REALSXP,1));
    REAL(res)[0] = logpkwmix<double>(Rf_asReal(x), Rf_asReal(logit_p), Rf_asReal(df), Rf_asReal(a), Rf_asReal(b));
    UNPROTECT(1);
    return res;
  }

  SEXP logdkwmixE0(SEXP x, SEXP logit_p, SEXP df, SEXP a, SEXP b){
    SEXP res;    
    PROTECT(res = Rf_allocVector(REALSXP,1));
    REAL(res)[0] = logdkwmixE0<double>(Rf_asReal(x), Rf_asReal(logit_p), Rf_asReal(df), Rf_asReal(a), Rf_asReal(b));
    UNPROTECT(1);
    return res;
  }
  SEXP logpkwmixE0(SEXP x, SEXP logit_p, SEXP df, SEXP a, SEXP b){
    SEXP res;
    PROTECT(res = Rf_allocVector(REALSXP,1));
    REAL(res)[0] = logpkwmixE0<double>(Rf_asReal(x), Rf_asReal(logit_p), Rf_asReal(df), Rf_asReal(a), Rf_asReal(b));
    UNPROTECT(1);
    return res;
  }
}
