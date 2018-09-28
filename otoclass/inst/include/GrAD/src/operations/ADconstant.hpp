// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADconstant : public ADnode<T>{

  T xval;
  
  ADconstant(T x0);
  
  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);
 
};


template<class T>
ADconstant<T>::ADconstant(T x0) :
  ADnode<T>("CONSTANT"), xval(x0) {
  this->setValue(x0);
}

template<class T>
T ADconstant<T>::fn(vector<T> x){
     return xval;
  }

template<class T>
vector<T> ADconstant<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  for(int i = 0; i < x.size(); ++i){
      res[i] = T(0.0);
  }
  return res;
}



template<class T>
void ADconstant<T>::bdfn(T w, vector<T>& theta){
  return;
}


template<class T>
void ADconstant<T>::update(vector<T>& theta){
  return;
}

