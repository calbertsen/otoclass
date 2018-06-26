#include "../inst/include/polygon.hpp"

extern "C" {
SEXP polygon_area(SEXP x, SEXP y){
  SEXP res;
  PROTECT(res = Rf_allocVector(REALSXP, 1));  
  double* pr = REAL(res);
  pr[0] = 0.0;
  double* px = REAL(x);
  double* py = REAL(y);
  // Shoelace formula
  for(int i = 1; i < Rf_length(x); ++i){
    //REAL(res)[0] += REAL(x)[i-1] * REAL(y)[i] - REAL(x)[i] * REAL(y)[i-1];
    pr[0] += px[i-1] * py[i] - px[i] * py[i-1];
  }
  REAL(res)[0] *= 0.5;
  UNPROTECT(1);
  return res;
}
}
