#include "convert.hpp"

extern "C" {
  SEXP biascorrect_logistic(SEXP x);
  SEXP biascorrect_logistic_gradient(SEXP x);
  SEXP biascorrect_objective(SEXP x, SEXP X, SEXP M);
  SEXP biascorrect_gradient(SEXP x, SEXP X, SEXP M);
  
}
