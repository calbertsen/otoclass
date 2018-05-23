#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include "convert.hpp"
// #include <Eigen/Dense>
using Eigen::Matrix;
// using Eigen::Dynamic;

// #include <Rinternals.h>

// Matrix<double, Dynamic, Dynamic> asDoubleMatrix(SEXP x);

template<class Type>
Matrix<Type, Dynamic, Dynamic> efd2coord_work(Matrix<Type, Dynamic, Dynamic> efd, int N);

extern "C" {
  SEXP efd2coordSEXP(SEXP efd, SEXP N);
}
