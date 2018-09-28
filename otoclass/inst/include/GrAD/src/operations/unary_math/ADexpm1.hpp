// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADexpm1 : public ADnode<T>{


  ADexpm1(ADnode<T>* L);
  ADexpm1(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADexpm1<T>::ADexpm1(ADnode<T>* L): ADnode<T>("ADEXPM1",L){
  T cv = expm1(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADexpm1<T>::ADexpm1(ADnode<T> L): ADnode<T>("ADEXPM1",&L){
  T cv = expm1(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADexpm1<T>::fn(vector<T> x){
  return expm1(this->ptrL->fn(x));
}

template<class T>
vector<T> ADexpm1<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * exp(fnL);
  }
  return res;      
}


template<class T>
AD<T> expm1(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADexpm1<T>(orx));
  return newAD;
}


template<class T>
void ADexpm1<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * exp(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADexpm1<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = expm1(this->ptrL->curval);
  this->setValue(cv);
  return;
}

