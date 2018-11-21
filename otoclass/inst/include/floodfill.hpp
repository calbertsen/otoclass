#include <vector>
#include "convert.hpp"

extern "C" {
  SEXP scanlineFill(SEXP pic, SEXP tol);
}
