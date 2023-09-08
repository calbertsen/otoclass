#include "../inst/include/polygon.hpp"

#define ASSERT(x,msg) if(!(x)) Rf_error(msg);

extern "C" {
SEXP polygon_area(SEXP x, SEXP y){
  double r = 0.0;
  double* px = REAL(x);
  double* py = REAL(y);
  // Shoelace formula
  for(int i = 1; i < Rf_length(x); ++i){
    //REAL(res)[0] += REAL(x)[i-1] * REAL(y)[i] - REAL(x)[i] * REAL(y)[i-1];
    r += px[i-1] * py[i] - px[i] * py[i-1];
  }
  r *= 0.5;
  return Rf_ScalarReal(r);
}

  // Shimrat (1962) Algorithm 112: Position of a point relative to polygon. DOI: 10.1145/368637.368653
  // Hacker (1962) Certification of Algorithm 112: Position of point relative to polygon. DOI: 10.1145/355580.369118 
  SEXP point_in_polygon(SEXP x0_, SEXP y0_, SEXP x, SEXP y){
    int n = Rf_length(x);
    double* px = REAL(x);
    double* py = REAL(y);
    double x0 = Rf_asReal(x0_);
    double y0 = Rf_asReal(y0_);
    bool b = true;
    for(int i = 0; i < n - 1; ++i){
      if( ((y0 <= py[i]) == (y0 > py[i+1])) &&
	  ((x0 - px[i] - (y0 - py[i]) * (px[i+1] - px[i]) / (py[i+1] - py[i])) < 0) ){
	b = !b;
      }
    }
    return Rf_ScalarLogical((int)!b);
  }
  
}
