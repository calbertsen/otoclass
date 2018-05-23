#include "../inst/include/efd.hpp"

template<class Type>
Matrix<Type, Dynamic, Dynamic> efd2coord_work(Matrix<Type, Dynamic, Dynamic> efd, int N){
  Matrix<Type, Dynamic, Dynamic> coord(2,N);
  coord.setZero();

  for(int tt = 0; tt < N; ++tt){
    Type t = (Type)tt / (Type)(N-1);
    for(int i = 0; i < efd.cols(); ++i){
      coord(0,tt) += efd(0,i) * cos(2.0 * Type(i+1) * M_PI * t);
      coord(0,tt) += efd(1,i) * sin(2.0 * Type(i+1) * M_PI * t);
      coord(1,tt) += efd(2,i) * cos(2.0 * Type(i+1) * M_PI * t);
      coord(1,tt) += efd(3,i) * sin(2.0 * Type(i+1) * M_PI * t);
    }
  }
  return coord;
}


extern "C" {
  SEXP efd2coordSEXP(SEXP efd, SEXP N){
    return asSEXP(efd2coord_work(asDoubleMatrix(efd), asInteger(N)));
  }
}

