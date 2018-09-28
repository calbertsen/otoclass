// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADprod : public ADnode<T>{


  ADprod(ADnode<T>* L, ADnode<T>* R);
  ADprod(ADnode<T> L, ADnode<T> R);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADprod<T>::ADprod(ADnode<T>* L, ADnode<T>* R): ADnode<T>("ADPROD",L,R){
  T cv = this->ptrL->curval * this->ptrR->curval;
  this->setValue(cv);
}
template<class T>
ADprod<T>::ADprod(ADnode<T> L, ADnode<T> R): ADnode<T>("ADPROD",&L,&R){
  T cv = this->ptrL->curval * this->ptrR->curval;
  this->setValue(cv);
}

template<class T>
T ADprod<T>::fn(vector<T> x){
  return this->ptrL->fn(x) * this->ptrR->fn(x);
}

template<class T>
vector<T> ADprod<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  vector<T> dfnR = this->ptrR->dfn(x);
  T fnL = this->ptrL->fn(x);
  T fnR = this->ptrR->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i]*fnR + dfnR[i]*fnL;
  }
  return res;      
}

template<class T>
void ADprod<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T fnR = this->ptrR->getValue();
  T wL = w * fnR;
  T wR = w * fnL;
  this->ptrL->bdfn(wL,theta);
  this->ptrR->bdfn(wR,theta);
  return;
}


template<class T>
void ADprod<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  this->ptrR->update(theta);
  T cv = this->ptrL->curval * this->ptrR->curval;
  this->setValue(cv);
  return;
}
