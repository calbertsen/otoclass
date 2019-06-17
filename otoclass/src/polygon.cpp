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

  SEXP polygon_centroid(SEXP x, SEXP y){
    SEXP res;
    PROTECT(res = Rf_allocVector(REALSXP, 2));  
    double* pr = REAL(res);
    pr[0] = 0.0; pr[1] = 0.0;
    double* px = REAL(x);
    double* py = REAL(y);
    SEXP A;
    PROTECT(A = polygon_area(x,y));
    for(int i = 1; i < Rf_length(x); ++i){
      double tt = (px[i-1] * py[i] - px[i] * py[i-1]);
      pr[0] += (px[i] + px[i-1]) * tt;
      pr[1] += (py[i] + py[i-1]) * tt;
    }
    REAL(res)[0] *= 1.0 / (6.0 * REAL(A)[0]);
    REAL(res)[1] *= 1.0 / (6.0 * REAL(A)[0]);
    UNPROTECT(2);
    return res;

  }
}

