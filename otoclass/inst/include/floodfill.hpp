#include <vector>
#include "convert.hpp"

extern "C" {
  SEXP scanlineFill(SEXP pic, SEXP tol, SEXP x0, SEXP y0, SEXP setInitCol);
}
