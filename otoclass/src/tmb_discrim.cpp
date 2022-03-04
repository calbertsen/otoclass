#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#define TMB_SAFEBOUNDS
#include "../inst/include/efd.hpp"


#define mkString Rf_mkString
#define mkChar Rf_mkChar
#define allocVector Rf_allocVector
#define ScalarInteger Rf_ScalarInteger
#define ScalarLogical Rf_ScalarLogical
#define isNull Rf_isNull
#define install Rf_install
#define findVar Rf_findVar
#define setAttrib Rf_setAttrib
#include <TMB.hpp>

#include "../inst/include/param_types.hpp"
#include "../inst/include/pnorm.hpp"
#include "../inst/include/spline.hpp"

using namespace density;
using Eigen::SparseMatrix;


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

// logdrobust can't handle NA, since x will be observation - prediction (i.e., include parameters)
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

// Logit function
double logit(double x){
  return log(x) - log((double)1.0 - x);
}

template<class Type>
AD<Type> logit(AD<Type> x){
  return log(x) - log((AD<Type>)1.0 - x);
}


template<class Type>
vector<Type> logit(vector<Type> x){
  vector<Type> r(x.size());
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = logit(x(i));
  return r;
}


template<class Type>
matrix<Type> logit(matrix<Type> x){
  matrix<Type> r(x.rows(), x.cols());
  r.setZero();
  for(int i = 0; i < r.rows(); ++i)
      for(int j = 0; j < r.cols(); ++j)
	r(i,j) = logit(x(i,j));
  return r;
}

double logistic(double x){
  return 1.0 / (1.0 + exp(-x));
}

template<class Type>
AD<Type> logistic(AD<Type> x){
  return 1.0 / (1.0 + exp(-x));
}
template<class Type>
vector<Type> logistic(vector<Type> x){
  vector<Type> r(x.size());
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = logistic(x(i));
  return r;
}



double log2logit(double logx){
  return logx - logspace_sub(0.0,logx);
}

template<class Type>
Type log2logit(Type logx){
  return logx - logspace_sub((Type)0.0,logx);
}


template<class Type>
AD<Type> log2logit(AD<Type> logx){
  return logx - logspace_sub((AD<Type>)0.0,logx);
}
template<class Type>
vector<Type> log2logit(vector<Type> x){
  vector<Type> r(x.size());
  r.setZero();
  for(int i = 0; i < r.size(); ++i)
    r(i) = log2logit(x(i));
  return r;
}




template<class Type>
vector<Type> toLogProportion(vector<Type> x){
  vector<Type> r(x.size()+1);
  r.setZero();
  Type lps = 0.0;
  for(int i = 0; i < x.size(); ++i){
    r(i) = x(i);
    lps = logspace_add2(lps,x(i));
  }
  r -= lps;
  return r;
}




template<class Type>
matrix<Type> naiveReshape(matrix<Type> x, int r, int c){
  if(x.size() != r * c)
    Rf_error("Wrong size in naiveReshape");
  matrix<Type> y(r,c);
  for(int i = 0; i < x.size(); ++i)
    y(i) = x(i);
  return y;
}

template<class Type>
matrix<Type> naiveReshapeRow(matrix<Type> x, int r){
  return naiveReshape(x, r, x.size() / r);
}

template<class Type>
matrix<Type> naiveReshapeCol(matrix<Type> x, int c){
  return naiveReshape(x, x.size() / c, c);
}


template<class Type>
Type dmultinom2(vector<Type> x, vector<Type> logp, int give_log=0){
  Type logres = 0.0;
  Type xsum = 0.0;
  for(int i = 0; i < x.size(); ++i){
    logres += -lgamma(x(i) + 1.0) + x(i) * logp(i);
    xsum += x(i);
  }
  logres += lgamma(xsum + 1.0);
  //-lgamma(alpha.sum()) + lgamma(alpha).sum() + xcontr.sum();
  if(give_log)return logres; else return exp(logres);
}


template<class Type>
Type ddirichletmultinom(vector<Type> x, vector<Type> logp, Type logAlphaScale, int give_log=0){
  vector<Type> alpha = exp(logp + logAlphaScale);
  Type logres = 0.0;
  Type alphsum = 0.0;
  Type xsum = 0.0;
  for(int i = 0; i < x.size(); ++i){
    logres += lgamma(x(i) + alpha(i)) - lgamma(x(i) + 1.0) - lgamma(alpha(i));
    alphsum += alpha(i);
    xsum += x(i);
  }
  logres += lgamma(xsum + 1.0) + lgamma(alphsum) - lgamma(xsum + alphsum);
  //-lgamma(alpha.sum()) + lgamma(alpha).sum() + xcontr.sum();
  if(give_log)return logres; else return exp(logres);
}



template <class Type>
struct OBSERVATION_DENSITY {
  virtual Type operator()(vector<Type> x, vector<Type> mu, vector<Type> scale) = 0;
  //virtual ~OBSERVATION_DENSITY()=0;
};

template<class Type>
struct SNP1 : OBSERVATION_DENSITY<Type> {
  Type operator()(vector<Type> x, vector<Type> mu, vector<Type> scale){
    Type r = 0.0;
    for(int i = 0; i < x.size(); ++i){
      Type logp = -logspace_add2((Type)0.0, -mu(i));
      Type log1mp = -logspace_add2((Type)0.0, mu(i));
      if(isNA(x(i))){
	r += 0.0;
      }else if(CppAD::Integer(x(i)) == 0){
	r += 2.0 * logp;
      }else if(CppAD::Integer(x(i)) == 1){
	r += log(2.0) + logp + log1mp;
      }else if(CppAD::Integer(x(i)) == 2){
	r += 2.0 * log1mp;
      }else{
        Rf_error("Error in SNP data");
      }
    }
    return -r;
  }
};

template<class Type>
struct SNP2 : OBSERVATION_DENSITY<Type> {
  int nLoci;
  int nAllele;
  
  SNP2() : nLoci(0), nAllele(0) {};
  SNP2(int nLoci_, int nAllele_) : nLoci(nLoci_), nAllele(nAllele_) {};
  
  // x: vectr of length #Allele x # Loci observations
  // mu: vector of length #(Allele-1) x #Loci; 
  // scale: vector of length #Loci; 
  Type operator()(vector<Type> x, vector<Type> mu, vector<Type> scale){
    matrix<Type> xx = asMatrix(x, nAllele, nLoci);
    matrix<Type> mm = asMatrix(mu, nLoci, nAllele-1).transpose();
    Type r = 0.0;
    for(int i = 0; i < xx.cols(); ++i){
      vector<Type> logp = toLogProportion((vector<Type>)mm.col(i));
      vector<Type> y = xx.col(i);
      if(sum(y) == 0){
	r += 0.0;
      }else if(R_finite(asDouble(scale(i)))){
	r += ddirichletmultinom(y,logp,log(scale(i)),true);
      }else if(!R_finite(asDouble(scale(i)))){
	r += dmultinom2(y,logp,true);
      }
    }
    return -r;
  }
};


template<class Type>
struct SPL_AR1 : OBSERVATION_DENSITY<Type> {
  vector<Type> xv;
  vector<Type> knots;
  
  SCALE_t<AR1_t<N01<Type> > > d0;
  
  SPL_AR1() : xv(0), knots(0) {};
  SPL_AR1(int N,
	  vector<Type> knots_,
	  Type sd,
	  Type phi) : xv(N), knots(knots_), d0(SCALE(AR1(phi), sd)) {
    for(int i = 0; i < N; ++i)
      xv(i) = Type(i) / Type(N-1);    
  };

  Type operator()(vector<Type> x, vector<Type> mu, vector<Type> scale){
    vector<Type> r(x.size());
    r.setZero();
    vector<Type> rs(x.size());
    rs.setZero();
    int i0 = 0;
    int i1 = x.size() - 1;
    for(int i = 0; i < r.size(); ++i){
      r(i) = x(i) - bcspline(xv(i), knots, (vector<Type>)mu.segment(1,mu.size()-1)) - mu(0);
      rs(i) = (bcspline(xv(i), knots, (vector<Type>)scale.segment(1,scale.size()-1)) - scale(0));
      if(i > 0 && !isNA(x(i)) && isNA(x(i-1))) // First non-NA
	i0 = i;
      if(i < x.size()-1 && isNA(x(i+1)) && !isNA(x(i))) // Last non-NA
	i1 = i;
    }
    if(i0 == 0 && isNA(x(0)))	// All is missing
      return 0.0;
    return d0((vector<Type>)r.segment(i0, i1 - i0 + 1) / exp((vector<Type>)rs.segment(i0, i1 - i0 + 1))) + sum((vector<Type>)rs.segment(i0, i1 - i0 + 1)); 
  }
};



template<class Type>
struct FS1_AR1 : OBSERVATION_DENSITY<Type> {
  vector<Type> xv;
  int Ne;
  matrix<Type> basis;
  
  SCALE_t<AR1_t<N01<Type> > > d0;

  // Handle missing!
  
  FS1_AR1() : xv(0), Ne(0), basis(0,0) {};
  FS1_AR1(int N,
	  int Ne_,
	  Type sd,
	  Type phi) : xv(N), Ne(Ne_), basis(N,2 * Ne_), d0(SCALE(AR1(phi), sd)) {
    xv.setZero();
    basis.setZero();
    for(int i = 0; i < N; ++i){
      xv(i) = Type(i) / Type(N-1);
      for(int j = 0; j < Ne_; ++j){
	basis(i,2 * j) = cos(Type(2.0 * M_PI) * Type(j) * xv(i));
	basis(i,2 * j + 1) = cos(Type(2.0 * M_PI) * Type(j) * xv(i));
      }
    }
  };

  Type operator()(vector<Type> x, vector<Type> mu, vector<Type> scale){
    vector<Type> r(x.size());
    if(isNA(x(0)))		// Assuming all are missing
      return(0.0);
    vector<Type> mv = basis * (vector<Type>)mu.segment(1,mu.size()-1);
    for(int i = 0; i < r.size(); ++i)
      r(i) = x(i) - mv(i) - mu(0);
    return d0(r); 
  }
};


template<class Type>
struct SPL_CATEGORICAL : OBSERVATION_DENSITY<Type> {
  int Ncat;
  vector<Type> knots;
  
  SPL_CATEGORICAL() : knots(0) {};
  SPL_CATEGORICAL(int N,
	  vector<Type> knots_) : Ncat(N), knots(knots_) {};

  Type operator()(vector<Type> x, vector<Type> mu, vector<Type> scale){
    // x.size() == 1
    if(isNA(x(0)))
      return(0.0);
    int xi = CppAD::Integer(x(0));
    vector<Type> vv(Ncat-1);
    vv.setZero();
    for(int i = 0; i < vv.size(); ++i)
      vv(i) = bcspline(Type(i), knots, (vector<Type>)mu.segment(1,mu.size()-1)) - mu(0);
    vector<Type> logP = toLogProportion(vv);
    return -logP(xi);
  }
};


// template <class Type>
// struct DistVec : vector<OBSERVATION_DENSITY<Type>* > {
//   ~DistVec(){
//     for(int i = 0; i < this->size(); ++i)
//       delete this->operator()(i);    
//   }
//   DistVec(int n) : vector<OBSERVATION_DENSITY<Type>* >(n) {};
//   DistVec() : vector<OBSERVATION_DENSITY<Type>* >() {};
  
// };

template <class Type>
struct MVMIX_t : OBSERVATION_DENSITY<Type> {
  Type halfLogDetS;         
  vector<Type> p1;                  /*fraction t3*/
  vector<Type> df;
  matrix<Type> Sigma;       
  vector<Type> sd;
  matrix<Type> L_Sigma;
  matrix<Type> inv_L_Sigma;
public:
  MVMIX_t(){}
  MVMIX_t(matrix<Type> Sigma_, vector<Type> p1_, vector<Type> df_) : OBSERVATION_DENSITY<Type>() {
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
    // vector<Type> tmp = (x - mu) / scale;
    vector<Type> z = inv_L_Sigma*x;
    return -sum(logdrobust(z,p1,df))+halfLogDetS;// + sum(log(scale));
  }
  Type operator()(vector<Type> x, vector<Type> mu, vector<Type> scale){
    vector<Type> keep(x.size());
    keep.setConstant(1.0);
    int hasNA = 0;
    for(int i = 0; i < x.size(); ++i)
      if(isNA(x(i))){
	x(i) = 0.0;
	keep(i) = 0.0;
	hasNA++;
      }
    if(hasNA == 0){
      vector<Type> tmp = (x - mu) / scale;
      Type v1 = this->operator()(tmp);
      Type v2 = sum(log(scale));
      return v1+v2;
    }
    if(hasNA == x.size()){
      return 0.0;
    }

    matrix<Type> S = Sigma;
    vector<Type> not_keep = Type(1.0) - keep;
    for(int i = 0; i < S.rows(); i++){
      for(int j = 0; j < S.cols(); j++){
  	S(i,j) = S(i,j) * keep(i) * keep(j);
      }
      //S(i,i) += not_keep(i) * pow((Type(1)-p1)*sqrt(Type(0.5)/M_PI)+p1*(Type(1)/M_PI),2); //(t(1))
      S(i,i) += not_keep(i) * pow((Type(1)-p1(i))*sqrt(Type(0.5)/M_PI)+p1(i)*(Type(2)/(M_PI*sqrt(df(i)))),2);
    }
    return MVMIX_t<Type>(S,p1,df)(x * keep, mu, scale);
  }

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
    MVMIX,
    SNP1,
    SNP2,
    SPL_AR1,
    FS1_AR1,
    SPL_CATEGORICAL
  };
}



// namespace Model{
//   enum Types {
//     Independent,
//     Unstructured,
//     GMRF,
//     RotatingOU
//   };
// }


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


template<class Type>
struct MatrixList : vector<matrix<Type> > {
  
  MatrixList() : vector<matrix<Type> >() {};
  MatrixList(int n) : vector<matrix<Type> >(n) {};
  MatrixList(SEXP x){ 
    (*this).resize(Rf_length(x));
    for(int i=0; i<Rf_length(x); i++){
      SEXP m = PROTECT(VECTOR_ELT(x, i));
      (*this)(i) = asMatrix<Type>(m);
      UNPROTECT(1);
    }
  }
};

template<class Type>
struct SparseMatrixList : vector<SparseMatrix<Type> > {
  
  SparseMatrixList() : vector<SparseMatrix<Type> >() {};
  SparseMatrixList(int n) : vector<SparseMatrix<Type> >(n) {};
  SparseMatrixList(SEXP x){ 
    (*this).resize(Rf_length(x));
    for(int i=0; i<Rf_length(x); i++){
      SEXP m = PROTECT(VECTOR_ELT(x, i));
      (*this)(i) = asSparseMatrix<Type>(m);
      UNPROTECT(1);
    }
  }
};


template <class Type>
Type objective_function<Type>::operator() () {
  DATA_INTEGER(model);
  Model::Types modelType = static_cast<Model::Types>(model);
    
  DATA_IMATRIX(G);
  // DATA_IMATRIX(Gxtra);
  DATA_IVECTOR(Gnlevels);
  DATA_STRUCT(Gconversion, MatrixList);
  DATA_ARRAY(Y);
  DATA_SPARSE_MATRIX(X);
  DATA_SPARSE_MATRIX(XCom);
  DATA_SPARSE_MATRIX(XLogScale);
  DATA_SPARSE_MATRIX(XTheta);
  DATA_STRUCT(Z, SparseMatrixList);
  DATA_STRUCT(ZCom, SparseMatrixList);
  DATA_STRUCT(ZTheta, SparseMatrixList);
  
  DATA_IVECTOR(penalty);

  DATA_ARRAY(Y_pred);
  DATA_IMATRIX(G_pred);
  DATA_SPARSE_MATRIX(X_pred);
  DATA_SPARSE_MATRIX(XCom_pred);
  DATA_SPARSE_MATRIX(XLogScale_pred);
  DATA_SPARSE_MATRIX(XTheta_pred);

  DATA_STRUCT(Z_pred, SparseMatrixList);
  DATA_STRUCT(ZCom_pred, SparseMatrixList);
  DATA_STRUCT(ZTheta_pred, SparseMatrixList);


  DATA_INTEGER(increaseFirstCoordinate);
  DATA_IMATRIX(identifyMatrix);

  DATA_VECTOR(knots);
  DATA_INTEGER(Ncat);
  
  PARAMETER_ARRAY(mu);
  PARAMETER_MATRIX(commonMu);
  PARAMETER_MATRIX(logSigma);
  PARAMETER_MATRIX(corpar);

  // For mu
  PARAMETER_CMOE_3DARRAY(U);
  PARAMETER_CMOE_3DARRAY(Ucor);
  PARAMETER_CMOE_3DARRAY(UlogSd);
  // For muCommon
  PARAMETER_CMOE_MATRIX(UCom);
  PARAMETER_CMOE_MATRIX(UComcor);
  PARAMETER_CMOE_MATRIX(UComlogSd);
  // For Theta
  PARAMETER_CMOE_MATRIX(UTheta);
  PARAMETER_CMOE_MATRIX(UThetacor);
  PARAMETER_CMOE_MATRIX(UThetalogSd);
  PARAMETER_VECTOR(logLambda);
  // PARAMETER_MATRIX(thetaIn);
  PARAMETER_ARRAY(betaLogScale);
  PARAMETER_MATRIX(betaTheta);
  PARAMETER_CMOE_MATRIX(MIn);
  PARAMETER_MATRIX(tmixpIn);
  PARAMETER_MATRIX(logDf);

  ////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// Check //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if(G.rows() != Y.cols())
    Rf_error("Wrong input length. G and Y must match");

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// Prepare /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  int nFeatures = 0;
  if(modelType == Model::MVMIX){
    nFeatures = Y.rows();
  }else if(modelType == Model::SNP1){
    nFeatures = Y.rows();
  }else if(modelType == Model::SNP2){
    nFeatures = betaLogScale.dim[1];
  }else if(modelType == Model::SPL_AR1){
    nFeatures = betaLogScale.dim[1];
  }else if(modelType == Model::FS1_AR1){
    nFeatures = betaLogScale.dim[1];
  }else if(modelType == Model::SPL_CATEGORICAL){
    nFeatures = betaLogScale.dim[1];
  }
  
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
	if(increaseFirstCoordinate && i == 0 && j == 0){
	  muUse(i,j,k) = muUse(i,j,k-1) + exp(mu(i,j,k));
	}else{
	  muUse(i,j,k) = muUse(i,j,0) + mu(i,j,k);
	}
      }
    }

  for(int i = 0; i < identifyMatrix.rows(); ++i){
    muUse(0,0,identifyMatrix(i,1)) = muUse(0,0,identifyMatrix(i,0)) + exp(muUse(0,0,identifyMatrix(i,1)));
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
    if(modelType != Model::SPL_AR1 &&
       modelType != Model::FS1_AR1 &&
       modelType != Model::SPL_CATEGORICAL){
      if(stmp.cols() > 1){
	stmp = UNSTRUCTURED_CORR((vector<Type>)corpar.col(k)).cov();
      }
      for(int i = 0; i < sigma.rows(); ++i)
	for(int j = 0; j < sigma.rows(); ++j)
	  stmp(i,j) *= sigma(i,k) * sigma(j,k);
      SigmaList(k) = stmp;
    }
  }

  
  // Transform theta
  // matrix<Type> theta(thetaIn.rows() + 1,thetaIn.cols());
  // for(int j = 0; j < theta.cols(); ++j){
  //   for(int i = 0; i < thetaIn.rows(); ++i)
  //     theta(i,j) = exp(thetaIn(i,j));
  //   theta(theta.rows()-1,j) = 1.0;
  //   theta.col(j) /= theta.col(j).sum();
  // }
  // matrix<Type> theta_logit = logit(theta);

  // vector<matrix<Type> > theta_reduced(Gconversion.size());
  // vector<matrix<Type> > theta_reduced_logit(Gconversion.size());
  // for(int i = 0; i < theta_reduced.size(); ++i){
  //   matrix<Type> tmp = Gconversion(i) * theta;
  //   theta_reduced(i) = tmp;
  //   theta_reduced_logit(i) = logit(tmp);
  // }
 

  // Transform confusionw
  vector<matrix<Type> > Mvec(MIn.cols());
  for(int k = 0; k < MIn.cols(); ++k){
    matrix<Type> Mtmp1 = MIn.col(k);
    matrix<Type> Mtmp(Mtmp1.rows(), Mtmp1.cols());
    Mtmp.setZero();
    for(int j = 0; j < Mtmp.cols(); ++j){
      for(int i = 0; i < Mtmp.rows(); ++i){
	if(i != j){
	  Mtmp(i,j) = exp(Mtmp1(i,j));
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
  std::vector<OBSERVATION_DENSITY<Type>* > dist(Gnlevels(0)); //(NLEVELS(G));
  for(int i = 0; i < Gnlevels(0); ++i){
    if(modelType == Model::MVMIX){
      dist[i] = new MVMIX_t<Type>(SigmaList(i), (vector<Type>)tmixp.col(i), (vector<Type>)df.col(i));
    }else if(modelType == Model::SNP1){
      dist[i] = new SNP1<Type>();
    }else if(modelType == Model::SNP2){
      dist[i] = new SNP2<Type>(Y.dim[1],Y.dim[0]);
    }else if(modelType == Model::SPL_AR1){
      dist[i] = new SPL_AR1<Type>(Y.dim[0],knots, exp(logSigma(0,i)), 1.0 / (1.0 + exp(-logSigma(1,i))));
      REPORT(knots);
    }else if(modelType == Model::FS1_AR1){
      dist[i] = new FS1_AR1<Type>(Y.dim[0],(nFeatures-1)/2, exp(logSigma(0,i)), 1.0 / (1.0 + exp(-logSigma(1,i))));
    }else if(modelType == Model::SPL_CATEGORICAL){
      dist[i] = new SPL_CATEGORICAL<Type>(Ncat,knots);
      REPORT(knots);
    }else{
      Rf_error("Model not implemented.");
    }
  }

  // Likelihood value
  Type nll = 0.0;

  // Random effects
  // U
  for(int i = 0; i < U.cols(); ++i){ // Loop over RE terms
    array<Type> Ui = U.col(i);
    for(int g = 0; g < Ui.cols(); ++g){ // Loop over groups
      for(int j = 0; j < Ui.col(g).cols(); ++j){ // Loop over features
	UNSTRUCTURED_CORR_t<Type> Uigjdens1(Ucor.col(i).col(g).col(j));
	vector<Type> UigjScale = exp((vector<Type>)UlogSd.col(i).col(g).col(j));
	int nr = UigjScale.size();
	VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > UigjDens = density::VECSCALE(Uigjdens1, UigjScale);
    	matrix<Type> Uigj = naiveReshapeRow((matrix<Type>)Ui.col(g).col(j), nr);	
	for(int k = 0; k < Uigj.cols(); ++k) // Loop over RE
	  nll += UigjDens((vector<Type>)Uigj.col(k));
      }
    }
  }

  // Ucommon
  for(int i = 0; i < UCom.cols(); ++i){ // Loop over RE terms
    matrix<Type> UCi = UCom.col(i);
    for(int j = 0; j < UCi.cols(); ++j){ // Loop over features
      UNSTRUCTURED_CORR_t<Type> UCijdens1(UComcor.col(i).col(j));
      vector<Type> UCijScale = exp((vector<Type>)UComlogSd.col(i).col(j));
      int nr = UCijScale.size();
      VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > UCijDens = density::VECSCALE(UCijdens1, UCijScale);
      matrix<Type> UCij = naiveReshapeRow((matrix<Type>)UCi.col(j),nr);
      for(int k = 0; k < UCij.cols(); ++k) // Loop over RE
  	nll += UCijDens((vector<Type>)UCij.col(k));
    }
  }


  // Utheta
  for(int i = 0; i < UTheta.cols(); ++i){ // Loop over RE terms
    matrix<Type> UTi = UTheta.col(i);
    for(int g = 0; g < UTi.cols(); ++g){ // Loop over groups
      UNSTRUCTURED_CORR_t<Type> UTigdens1(UThetacor.col(i).col(g));
      vector<Type> UTigScale = exp((vector<Type>)UThetalogSd.col(i).col(g));
      VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > UTigDens = density::VECSCALE(UTigdens1, UTigScale);
      int nr = UTigScale.size();
      matrix<Type> UTig = naiveReshapeRow((matrix<Type>)UTi.col(g),nr);
      for(int j = 0; j < UTig.cols(); ++j) // Loop over RE 	
  	nll += UTigDens((vector<Type>)UTig.col(j));
    }
  }


  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Estimate /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  array<Type> posterior_mean(nFeatures, Y.cols(), Gnlevels(0));
  posterior_mean.setZero();
  matrix<Type> prior_logitprobability(Gnlevels(0), Y.cols());
  prior_logitprobability.setZero();
  matrix<Type> posterior_logprobability_shape(Gnlevels(0), Y.cols());
  posterior_logprobability_shape.setZero();
  matrix<Type> posterior_logprobability_all(Gnlevels(0), Y.cols());
  posterior_logprobability_all.setZero();
  matrix<Type> posterior_logprobability_np(Gnlevels(0), Y.cols());
  posterior_logprobability_np.setZero();
  for(int i = 0; i < Y.cols(); ++i){ // Loop over individuals
    //vector<Type> th = thetaDisp.col(dispersionGroup(i)).col(proportionGroup(i));
    vector<Type> thtmp(Gnlevels(0)-1);
    thtmp.setZero();

    for(int j = 0; j < thtmp.size(); ++j){
      thtmp(j) = ((vector<Type>)XTheta.col(i) * (vector<Type>)betaTheta.col(j)).sum();
    }
    for(int zz = 0; zz < ZTheta.size(); ++zz){
      thtmp += (vector<Type>)(ZTheta(zz).col(i).transpose() *  UTheta.col(zz));
    }
    vector<Type> logTh = toLogProportion(thtmp);
    // if(isNA(dispersionGroup(i))){
    //   th = theta.col(proportionGroup(i));
    // }else{
    //   for(int qq = 0; qq < thetaIn.rows(); ++qq)
    // 	th(qq) = exp(thetaIn(qq,proportionGroup(i)) + dispersion(qq,dispersionGroup(i)));
    //   th(th.size() - 1) = 1.0;
    //   th /= th.sum();
    // }
    prior_logitprobability.col(i) = log2logit(logTh);
    // vector<Type> Muse = Mvec(confusionGroup(i)).row(G(i));
    // Type NormalizationConstant = (Muse*th).sum();
    Type lps = R_NegInf;
    Type postProbNorm = R_NegInf;
    Type postProbNormNP = R_NegInf;
    // Contribution from G
    // nll -= log(NormalizationConstant);
    // Contrbution from Y
    for(int j = 0; j < Gnlevels(0); ++j){ // Loop over groups
      vector<Type> tmp = Y.col(i);
      vector<Type> tmpMean = (vector<Type>)(X.col(i).transpose() * muUse.col(j).matrix());
      for(int zz = 0; zz < Z.size(); ++zz){
	tmpMean += (vector<Type>)(Z(zz).col(i).transpose() * U.col(zz).col(j).matrix());
      }
      if(XCom.rows() > 0){
	tmpMean += (vector<Type>)(XCom.col(i).transpose() * commonMu);
      }
      for(int zz = 0; zz < ZCom.size(); ++zz){
	tmpMean += (vector<Type>)(ZCom(zz).col(i).transpose() * UCom.col(zz));
      }
      vector<Type> scale(nFeatures);
      scale.setConstant(1.0);
      if(modelType == Model::MVMIX || modelType == Model::SNP2){
	scale = exp((vector<Type>)(XLogScale.col(i).transpose() * betaLogScale.col(j).matrix()));
      }else if(modelType == Model::SPL_AR1){
	scale = (vector<Type>)(XLogScale.col(i).transpose() * betaLogScale.col(j).matrix());
      }
      // tmp -= tmpMean;
      if(modelType != Model::SNP2)
	posterior_mean.col(j).col(i) = tmpMean;
      // P(Y | S)
      Type ld = -dist[j]->operator()(tmp,tmpMean,scale); // dist returns negative log-likelihood!
      // P(Y | S) * P(S)      
      Type ld2 = ld + logTh(j);

      // P(C_i | S)
      for(int cc = 0; cc < G.cols(); ++cc){
	if(!isNA(G(i,cc))){
	  matrix<Type> MvecUse = Mvec(cc) * Gconversion(cc); // Columns should sum to one
	  ld2 += log(MvecUse(G(i,cc),j));
	}
      }

      // lp(j) = ld + log(th(j)) + log(Muse(j)) - log(NormalizationConstant);
      posterior_logprobability_shape(j,i) = ld + logTh(j);
      posterior_logprobability_all(j,i) = ld2;
      posterior_logprobability_np(j,i) = ld;
      lps = logspace_add2(lps,ld2);
      postProbNorm = logspace_add2(postProbNorm, posterior_logprobability_shape(j,i));
      postProbNormNP = logspace_add2(postProbNormNP, ld);
    }
    for(int j = 0; j < Gnlevels(0); ++j){
      posterior_logprobability_shape(j,i) -= postProbNorm;
      posterior_logprobability_all(j,i) -= lps;
      posterior_logprobability_np(j,i) -= postProbNormNP;
    }
    nll -= lps;
  }

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// Regularization ///////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if(penalty(0) != 0){
    for(int i = 0; i < mu.dim[0]; ++i)
      for(int j = 0; j < mu.dim[1]; ++j){
	for(int k = 0; k < mu.dim[2]; ++k){

	  // No penalty on intercept for first group
	  if(k == 0 && i == 0)
	    continue;
	  
	  Type sdp = 1.0;
	  if(modelType == Model::MVMIX)
	    sdp = sigma(j,k);
	  
	  switch(penalty(0)){
	  case -1:
	    nll -= dt(mu(i,j,k) / lambda(0),Type(3.0),true) - log(lambda(0));
	    break;
	  case 1:
	    nll -= -fabs(mu(i,j,k)) / sdp / lambda(0) - log(2.0 * lambda(0));
	    break;
	  case 2:	    
	    nll -= dnorm(mu(i,j,k),Type(0.0),lambda(0) * sdp, true);
	    break;
	  default:
	    nll -= my_atomic::Lp(mu(i,j,k),(Type)penalty(0),lambda(0) * sdp);
	  }
	}
      }
  }
  if(penalty(1) != 0){
    vector<Type> sigmaMeans(commonMu.cols());
    sigmaMeans.setConstant(1.0);
    if(modelType == Model::MVMIX)
      sigmaMeans = sigma.rowwise().mean();
    
    for(int i = 0; i < commonMu.rows(); ++i)
      for(int j = 0; j < commonMu.cols(); ++j){
	switch(penalty(1)){
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
	  nll -= my_atomic::Lp(commonMu(i,j),(Type)penalty(1),lambda(1) * sigmaMeans(j));
	}
      }
  }

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// Predict //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  if(Y_pred.cols() == 0){	// Group mean/proportion
    if(X_pred.cols() > 0){	// Group mean
      array<Type> GROUP_mean_group(muUse.col(0).cols(), X_pred.cols(), Gnlevels(0));
      GROUP_mean_group.setZero();
      matrix<Type> GROUP_mean_common(commonMu.cols(), X_pred.cols());
      GROUP_mean_common.setZero();
      array<Type> GROUP_mean_total(muUse.col(0).cols(), X_pred.cols(), Gnlevels(0));
      GROUP_mean_total.setZero();
      
      for(int i = 0; i < X_pred.cols(); ++i){
	if(XCom_pred.rows() > 0 && XCom_pred.cols() > 0)
	  GROUP_mean_common.col(i) = (vector<Type>)(XCom_pred.col(i).transpose() * commonMu);
	
	for(int j = 0; j < Gnlevels(0); ++j){
	  GROUP_mean_group.col(j).col(i) = (vector<Type>)(X_pred.col(i).transpose() * muUse.col(j).matrix());
	  GROUP_mean_total.col(j).col(i) =  (vector<Type>)GROUP_mean_group.col(j).col(i) + (vector<Type>)GROUP_mean_common.col(i);
	}
      }
      REPORT(GROUP_mean_group);
      REPORT(GROUP_mean_common);
      REPORT(GROUP_mean_total);
      ADREPORT(GROUP_mean_group);
      ADREPORT(GROUP_mean_common);
      ADREPORT(GROUP_mean_total);
    }
    if(XTheta_pred.cols() > 0){	// Group proportion
      matrix<Type> GROUP_logitprobability(Gnlevels(0), XTheta_pred.cols());
      GROUP_logitprobability.setZero();
      for(int i = 0; i < XTheta_pred.cols(); ++i){
	vector<Type> thtmp(Gnlevels(0)-1);
	thtmp.setZero();
	for(int j = 0; j < thtmp.size(); ++j){
	  thtmp(j) = ((vector<Type>)XTheta_pred.col(i) * (vector<Type>)betaTheta.col(j)).sum();
	}
	for(int zz = 0; zz < ZTheta_pred.size(); ++zz){
	  thtmp += (vector<Type>)(ZTheta_pred(zz).col(i).transpose() *  UTheta.col(zz));
	}  
	vector<Type> logTh = toLogProportion(thtmp);
	GROUP_logitprobability.col(i) = log2logit(logTh);
      }
      REPORT(GROUP_logitprobability);
      ADREPORT(GROUP_logitprobability);
    }
    
  }else{
    array<Type> pred_posterior_mean(nFeatures, Y_pred.cols(), Gnlevels(0));
    pred_posterior_mean.setZero();
    matrix<Type> pred_prior_logitprobability(Gnlevels(0), Y_pred.cols());
    pred_prior_logitprobability.setZero();
    matrix<Type> pred_posterior_logprobability_shape(Gnlevels(0), Y_pred.cols());
    pred_posterior_logprobability_shape.setZero();
    matrix<Type> pred_posterior_logprobability_all(Gnlevels(0), Y_pred.cols());
    pred_posterior_logprobability_all.setZero();
    matrix<Type> pred_posterior_logprobability_np(Gnlevels(0), Y_pred.cols());
    pred_posterior_logprobability_np.setZero();
  
    for(int i = 0; i < Y_pred.cols(); ++i){ // Loop over individuals
      vector<Type> thtmp(Gnlevels(0)-1);
      thtmp.setZero();
      
      for(int j = 0; j < thtmp.size(); ++j){
	thtmp(j) = ((vector<Type>)XTheta_pred.col(i) * (vector<Type>)betaTheta.col(j)).sum();
      }
      for(int zz = 0; zz < ZTheta_pred.size(); ++zz){
	thtmp += (vector<Type>)(ZTheta_pred(zz).col(i).transpose() *  UTheta.col(zz));
      }
      vector<Type> logTh = toLogProportion(thtmp);
      pred_prior_logitprobability.col(i) = log2logit(logTh);
    
      Type lps = R_NegInf;
      Type postProbNorm = R_NegInf;
      Type postProbNormNP = R_NegInf;
  
      for(int j = 0; j < Gnlevels(0); ++j){ // Loop over groups
	vector<Type> tmp = Y_pred.col(i);
	vector<Type> tmpMean = (vector<Type>)(X_pred.col(i).transpose() * muUse.col(j).matrix());
	for(int zz = 0; zz < Z_pred.size(); ++zz){
	  tmpMean += (vector<Type>)(Z_pred(zz).col(i).transpose() * U.col(zz).col(j).matrix());
	}
	if(XCom_pred.rows() > 0){
	  tmpMean += (vector<Type>)(XCom_pred.col(i).transpose() * commonMu);
	}
	for(int zz = 0; zz < ZCom_pred.size(); ++zz){
	  tmpMean += (vector<Type>)(ZCom_pred(zz).col(i).transpose() * UCom.col(zz));
	}
	vector<Type> scale(nFeatures);
	scale.setConstant(1.0);
	if(modelType == Model::MVMIX || modelType == Model::SNP2){
	  scale = exp((vector<Type>)(XLogScale_pred.col(i).transpose() * betaLogScale.col(j).matrix()));
	}else if(modelType == Model::SPL_AR1){
	  scale = (vector<Type>)(XLogScale_pred.col(i).transpose() * betaLogScale.col(j).matrix());
	}
	// tmp -= tmpMean;
	if(modelType != Model::SNP2)
	  pred_posterior_mean.col(j).col(i) = tmpMean;
	// P(Y | S)
	Type ld = -dist[j]->operator()(tmp,tmpMean,scale); // dist returns negative log-likelihood!
	// P(Y | S) * P(S)      
	Type ld2 = ld + logTh(j);

	// P(C_i | S)
	for(int cc = 0; cc < G_pred.cols(); ++cc){
	  if(!isNA(G_pred(i,cc))){
	    matrix<Type> MvecUse = Mvec(cc) * Gconversion(cc); // Columns should sum to one
	    ld2 += log(MvecUse(G_pred(i,cc),j));
	  }
	}
      
	pred_posterior_logprobability_shape(j,i) = ld + logTh(j);
	pred_posterior_logprobability_all(j,i) = ld2;
	pred_posterior_logprobability_np(j,i) = ld;
	lps = logspace_add2(lps,ld2);
        postProbNorm = logspace_add2(postProbNorm, pred_posterior_logprobability_shape(j,i));
	postProbNormNP = logspace_add2(postProbNormNP, ld);
      }
      for(int j = 0; j < Gnlevels(0); ++j){
	pred_posterior_logprobability_shape(j,i) -= postProbNorm;
	pred_posterior_logprobability_all(j,i) -= lps;
	pred_posterior_logprobability_np(j,i) -= postProbNormNP;      
      }
    }
    REPORT(pred_posterior_mean);
    REPORT(pred_posterior_logprobability_shape);
    REPORT(pred_posterior_logprobability_all);
    REPORT(pred_posterior_logprobability_np);
    REPORT(pred_prior_logitprobability);
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// Output //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  REPORT(muUse);
  ADREPORT(muUse);

  REPORT(commonMu);
  ADREPORT(commonMu);
    
  REPORT(posterior_mean);
  REPORT(posterior_logprobability_shape);
  REPORT(posterior_logprobability_all);
  REPORT(posterior_logprobability_np);
  REPORT(prior_logitprobability);

  // ADREPORT(posterior_mean);
  // ADREPORT(posterior_logprobability_shape);
  // ADREPORT(posterior_logprobability_all);


  // ADREPORT(pred_posterior_mean);
  // ADREPORT(pred_posterior_logprobability_shape);
  // ADREPORT(pred_posterior_logprobability_all);


  // REPORT(theta);
  // ADREPORT(theta);

  // REPORT(theta_logit);
  // ADREPORT(theta_logit);

  // REPORT(theta_reduced);
  // REPORT(theta_reduced_logit);
  
  REPORT(Mvec);

  REPORT(SigmaList);

  REPORT(tmixp);

  REPORT(df);


  
  return nll;
}
