#include "../inst/include/convert.hpp"

// From TMB: https://github.com/kaskr/adcomp/blob/master/TMB/inst/include/atomic_convolve.hpp
Eigen::MatrixXd convol2d_work(const Eigen::MatrixXd& x,
                              const Eigen::MatrixXd& K) {
  int kr = K.rows();
  int kc = K.cols();
  Eigen::MatrixXd y(x.rows() - kr + 1, x.cols() - kc + 1);
  for (int i = 0; i < y.rows(); i++)
    for (int j = 0; j < y.cols(); j++)
      y(i, j) = (x.block(i, j, kr, kc).array() * K.array()).sum();
  return y;
}


extern "C" {
  SEXP convol2d(SEXP x, SEXP k){
    return asSEXP(convol2d_work(asDoubleMatrix(x), asDoubleMatrix(k)));
  }

}
