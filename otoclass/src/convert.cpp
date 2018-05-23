// Modified from github.com/calbertsen/covafillr
#include "../inst/include/convert.hpp"

// Implemented in TMB
// SEXP asSEXP(const double &x){
//   SEXP y = PROTECT(Rf_allocVector(REALSXP,1));
//   REAL(y)[0] = x;

//   UNPROTECT(1);
//   return y; 

// }

// Implemented in TMB
// SEXP asSEXP(const int &x){
//   SEXP y = PROTECT(Rf_allocVector(INTSXP,1));
//   INTEGER(y)[0] = x;

//   UNPROTECT(1);
//   return y; 
// }

SEXP asSEXP(doubleVector x){
  int n = x.size();
  SEXP y = PROTECT(Rf_allocVector(REALSXP,n));

  for(int i = 0; i < n; ++i)
    REAL(y)[i] = x[i];

  UNPROTECT(1);
  return y; 
}

SEXP asSEXP(MatrixXd x){
  int nr = x.rows();
  int nc = x.cols();
  SEXP y = PROTECT(Rf_allocMatrix(REALSXP,nr,nc));
  for(int i = 0; i < nr; ++i)
    for(int j = 0; j < nc; ++j)
      REAL(y)[i+nr*j] = x(i,j);

  UNPROTECT(1);
  return y;
}

double asDouble(SEXP x){
  if(!(Rf_isNumeric(x) && Rf_length(x)==1))
    Rf_error("Element must be a numeric of length 1");
  return REAL(x)[0];
}

int asInteger(SEXP x){
  if(!(Rf_isInteger(x) && Rf_length(x)==1))
    Rf_error("Element must be an integer vector of length 1");
  return INTEGER(x)[0];
}

doubleVector asDoubleVector(SEXP x) {

  if(!Rf_isNumeric(x))
    Rf_error("Element must be a numeric vector");

  int n = Rf_length(x);
  doubleVector y(n);

  for(int i = 0; i < n; ++i)
    y(i) = REAL(x)[i];

  return y;
}


intVector asIntVector(SEXP x) {

  if(!Rf_isInteger(x))
    Rf_error("Element must be an integer vector");

  int n = Rf_length(x);
  intVector y(n);

  for(int i = 0; i < n; ++i)
    y(i) = INTEGER(x)[i];

  return y;
}


MatrixXd asDoubleMatrix(SEXP x)
{

  if(!Rf_isMatrix(x))
    Rf_error("Element must be a matrix");

  int nr = Rf_nrows(x);
  int nc = Rf_ncols(x);
  MatrixXd y(nr,nc);

  for(int i = 0; i < nr; ++i)
    for(int j = 0; j < nc; ++j)
      y(i,j)=REAL(x)[i+nr*j];

  return y;
}

