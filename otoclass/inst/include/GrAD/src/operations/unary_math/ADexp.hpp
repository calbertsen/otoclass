// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADexp : public ADnode<T>{


  ADexp(ADnode<T>* L);
  ADexp(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADexp<T>::ADexp(ADnode<T>* L): ADnode<T>("ADEXP",L){
  T cv = exp(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADexp<T>::ADexp(ADnode<T> L): ADnode<T>("ADEXP",&L){
  T cv = exp(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADexp<T>::fn(vector<T> x){
  return exp(this->ptrL->fn(x));
}

template<class T>
vector<T> ADexp<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * exp(fnL);
  }
  return res;      
}


template<class T>
AD<T> exp(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADexp<T>(orx));
  return newAD;
}

template<class T>
void ADexp<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * exp(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}

template<class T>
void ADexp<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = exp(this->ptrL->curval);
  this->setValue(cv);
  return;
}

