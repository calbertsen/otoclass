
#include "convert.hpp"

// From TMB: https://github.com/kaskr/adcomp/blob/master/TMB/inst/include/atomic_convolve.hpp
Eigen::MatrixXd convol2d_work(const Eigen::MatrixXd& x,
                              const Eigen::MatrixXd& K);

extern "C" {
  SEXP convol2d(SEXP x, SEXP k);
}
