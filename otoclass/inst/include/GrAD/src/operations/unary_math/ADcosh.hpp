// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADcosh : public ADnode<T>{


  ADcosh(ADnode<T>* L);
  ADcosh(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADcosh<T>::ADcosh(ADnode<T>* L): ADnode<T>("ADCOSH",L){
  T cv = cosh(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADcosh<T>::ADcosh(ADnode<T> L): ADnode<T>("ADCOSH",&L){
  T cv = cosh(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADcosh<T>::fn(vector<T> x){
  return cosh(this->ptrL->fn(x));
}

template<class T>
vector<T> ADcosh<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * sinh(fnL);
  }
  return res;      
}


template<class T>
AD<T> cosh(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADcosh<T>(orx));
  return newAD;
}


template<class T>
void ADcosh<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * sinh(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADcosh<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = cosh(this->ptrL->curval);
  this->setValue(cv);
  return;
}

