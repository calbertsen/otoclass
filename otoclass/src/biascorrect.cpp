// #include <iostream>

// template<class Type>
// class isDouble;

// class Rcout: public std::ostream;




#include "../inst/include/biascorrect.hpp"
#define NOUSEMEM
#include "../inst/include/GrAD/GrAD"



namespace biascorrect {

  using Eigen::Matrix;
  
  template<class Float>
  Array<Float,Dynamic,1> logistic(Array<Float,Dynamic,1> x){
    Array<Float,Dynamic,1> y(x.size()+1);
    Float ysum = 1.0;
    for(int i = 0; i < (int)x.size(); ++i){
      y(i) = exp(x(i));
      ysum += y(i);
    }
    y(y.size()-1) = 1.0;
    return y / ysum;
  }

  MatrixXd logistic_gr(doubleVector x){
    GrAD::ADparlist<double>* grd = new GrAD::ADparlist<double>();
    typedef GrAD::AD<double> AD;
    Array<AD,Dynamic,1> xf(x.size());
    for(int i = 0; i < (int)x.size(); ++i){
      xf(i) = AD(x(i));
      grd->Independent(xf(i));
    }
    
    Array<AD,Dynamic,1> y = logistic(xf);
    MatrixXd m(y.size(),x.size());
    m.setZero();    
    for(int i = 0; i < (int)y.size(); ++i){
      GrAD::vector<double> tmp = y(i).gr();//.getDeriv();
      for(int j = 0; j < (int)tmp.size(); ++j)
      	m(i,j) = tmp[j];
    }
    delete grd;
    return m;
  }

  
  template<class Float>
  Float fn(Array<Float,Dynamic,1> x, doubleVector X, MatrixXd M){
    Matrix<Float,Dynamic,Dynamic> pest(x.size()+1,1);
      pest = logistic(x);
      Matrix<Float,Dynamic,Dynamic> Mf = M.template cast<Float>();
    Array<Float,Dynamic,1> Xf = X.template cast<Float>();
    Array<Float,Dynamic,1> puse(Mf.rows());
    puse.setZero();
    for(int i = 0; i < Mf.cols(); ++i)
      for(int j = 0; j < Mf.rows(); ++j)
	puse(i) += Mf(i,j) * pest(j);
    Float nll = 0.0;
    for(int i = 0; i < puse.size(); ++i)
      nll -= Xf(i) * log(puse(i));
    return nll;
  
  }

  doubleVector gr(doubleVector x, doubleVector X, MatrixXd M){
  GrAD::ADparlist<double>* grd = new GrAD::ADparlist<double>();
  typedef GrAD::AD<double> AD;
    Array<AD,Dynamic,1> xf(x.size());
    for(int i = 0; i < (int)x.size(); ++i){
      xf(i) = AD(x(i));
      grd->Independent(xf(i));
    }
    
    AD nll = fn(xf,X,M);
    doubleVector r(x.size());
    r.setZero();    
    GrAD::vector<double> tmp = nll.gr();//.getDeriv();
    for(int j = 0; j < (int)tmp.size(); ++j)
      r(j) = tmp[j];
    delete grd;
    return r;
  }

}


extern "C" {

  SEXP biascorrect_logistic(SEXP x){
    return asSEXP(biascorrect::logistic(asDoubleVector(x)));
  }

  SEXP biascorrect_logistic_gradient(SEXP x){
    return asSEXP(biascorrect::logistic_gr(asDoubleVector(x)));
  }

  SEXP biascorrect_objective(SEXP x, SEXP X, SEXP M){
    return asSEXP(biascorrect::fn(asDoubleVector(x), asDoubleVector(X), asDoubleMatrix(M)));
  }

  SEXP biascorrect_gradient(SEXP x, SEXP X, SEXP M){
    return asSEXP(biascorrect::gr(asDoubleVector(x), asDoubleVector(X), asDoubleMatrix(M)));
  }

}
