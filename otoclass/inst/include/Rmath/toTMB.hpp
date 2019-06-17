#include "define.hpp"
#include "lbeta.hpp"
#include "pnorm.hpp"
#include "pt.hpp"


template<class Type>
Type pt(Type x, Type n, int lower_tail, int log_p){
  CppAD::vector<Type> tx(5);
  tx[0] = x;
  tx[1] = n;
  tx[2] = lower_tail;
  tx[3] = log_p;
  tx[4] = 0; // extra argument for derivative order
  Type res = my_atomic::pt2(tx)[0];
  return res;
}

#define VECTORIZE4_ttii(FUN)			\
  GVECTORIZE(FUN,V,T,I,I,N,N)			\
  GVECTORIZE(FUN,T,V,I,I,N,N)			\
  GVECTORIZE(FUN,V,V,I,I,N,N)


VECTORIZE4_ttii(pt)


template<class Type>
Type pnorm5(Type x, Type mu, Type sigma, int lower_tail, int log_p){
  CppAD::vector<Type> tx(6);
  tx[0] = x;
  tx[1] = mu;
  tx[2] = sigma;
  tx[3] = lower_tail;
  tx[4] = log_p;
  tx[5] = 0; // extra argument for derivative order
  Type res = my_atomic::pnorm5_2(tx)[0];
  return res;
}


#define VECTORIZE5_tttii(FUN)			\
  GVECTORIZE(FUN,V,T,T,I,I,N)			\
  GVECTORIZE(FUN,T,V,T,I,I,N)			\
  GVECTORIZE(FUN,T,T,V,I,I,N)			\
  GVECTORIZE(FUN,V,V,T,I,I,N)			\
  GVECTORIZE(FUN,T,V,V,I,I,N)			\
  GVECTORIZE(FUN,V,T,V,I,I,N)			\
  GVECTORIZE(FUN,V,V,V,I,I,N)

VECTORIZE5_tttii(pnorm5)
