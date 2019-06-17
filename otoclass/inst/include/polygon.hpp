#include "convert.hpp"

extern "C" {
  SEXP polygon_area(SEXP x, SEXP y);
  SEXP polygon_centroid(SEXP x, SEXP y);
}
