
// #define WITH_LIBTMB
// #include <TMB.hpp>
//#include "tmb_declarations_needed.hpp"

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <Eigen/Dense>

#define R_NO_REMAP
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::MatrixXd;

typedef Array<double,Dynamic, 1> doubleVector;
typedef Array<int, Dynamic, 1> intVector;

SEXP asSEXP(const double &x);
SEXP asSEXP(const int &x);
SEXP asSEXP(doubleVector x);
SEXP asSEXP(MatrixXd x);
double asDouble(SEXP x);
int asInteger(SEXP x);
doubleVector asDoubleVector(SEXP x);
intVector asIntVector(SEXP x);
MatrixXd asDoubleMatrix(SEXP x);
