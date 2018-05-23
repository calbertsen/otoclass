#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#include <Eigen/Dense>

template<class Type>
struct isDouble;

namespace CppAD {
  
template <class Type>
class vector;

}

using namespace Eigen;

template <class Type>
struct vector;

// {
//   typedef Type value_type;
//   typedef Array<Type,Dynamic,1> Base;
//   vector(void);

//   template<class T1>
//   vector(T1 x);

//   template<class T1, class T2>
//   vector(T1 x, T2 y);
  
//   template<class T1>
//   vector & operator= (const T1 & other);
  
//   using Base::operator();

//   vector<Type> operator()(vector<int> ind);

//   template<class T>
//   operator CppAD::vector<T>();

//   template<class T>
//   vector(CppAD::vector<T> x);

//   vector<Type> vec();
// };



template <class Type>
struct matrix;

// {
//   typedef Matrix<Type,Dynamic,Dynamic> Base;
//   matrix(void);
//   template<class T1>
//   matrix(T1 x);
//   template<class T1, class T2>
//   matrix(T1 x, T2 y);
//   template<class T1>
//   matrix & operator= (const T1 & other);
//   vector<Type> vec();
// };
