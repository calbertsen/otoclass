// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADusubtract : public ADnode<T>{


  ADusubtract(ADnode<T>* L);
  ADusubtract(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADusubtract<T>::ADusubtract(ADnode<T>* L): ADnode<T>("ADUSUBTRACT",L){
  T cv = -this->ptrL->curval;
  this->setValue(cv);
}
template<class T>
ADusubtract<T>::ADusubtract(ADnode<T> L): ADnode<T>("ADSUBTRACT",&L){
  T cv = -this->ptrL->curval;
  this->setValue(cv);
}

template<class T>
T ADusubtract<T>::fn(vector<T> x){
  T fnL = this->ptrL->fn(x);
  return -fnL;
}

template<class T>
vector<T> ADusubtract<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = -dfnL[i];
  }
  return res;      
}



template<class T>
void ADusubtract<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = -w;
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADusubtract<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = - this->ptrL->curval;
  this->setValue(cv);
  return;
}

