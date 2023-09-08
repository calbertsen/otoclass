#include "convert.hpp"



extern "C" {
  SEXP polygon_area(SEXP x, SEXP y);
  SEXP point_in_polygon(SEXP x0_, SEXP y0_, SEXP x, SEXP y);
}
