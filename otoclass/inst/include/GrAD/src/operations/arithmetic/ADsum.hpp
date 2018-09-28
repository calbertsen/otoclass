// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADsum : public ADnode<T>{


  ADsum(ADnode<T>* L, ADnode<T>* R);
  ADsum(ADnode<T> L, ADnode<T> R);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADsum<T>::ADsum(ADnode<T>* L, ADnode<T>* R): ADnode<T>("ADSUM",L,R){
  T cv = this->ptrL->curval + this->ptrR->curval;
  this->setValue(cv);
}
template<class T>
ADsum<T>::ADsum(ADnode<T> L, ADnode<T> R): ADnode<T>("ADSUM",&L,&R){
  T cv = this->ptrL->curval + this->ptrR->curval;
  this->setValue(cv);
}

template<class T>
T ADsum<T>::fn(vector<T> x){
  T fnL = this->ptrL->fn(x);
  T fnR = this->ptrR->fn(x);
  return fnL+fnR;
}

template<class T>
vector<T> ADsum<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  vector<T> dfnR = this->ptrR->dfn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] + dfnR[i];
  }
  return res;      
}


template<class T>
void ADsum<T>::bdfn(T w, vector<T>& theta){
  // T fnL = this->ptrL->getValue();
  // T fnR = this->ptrR->getValue();
  T wL = w;
  T wR = w;
  this->ptrL->bdfn(wL,theta);
  this->ptrR->bdfn(wR,theta);
  return;
}


template<class T>
void ADsum<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  this->ptrR->update(theta);
  T cv = this->ptrL->curval + this->ptrR->curval;
  this->setValue(cv);
  return;
}

