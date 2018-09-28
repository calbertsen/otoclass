// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADpow : public ADnode<T>{


  ADpow(ADnode<T>* L, ADnode<T>* R);
  ADpow(ADnode<T> L,ADnode<T> R);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);


};

template<class T>
ADpow<T>::ADpow(ADnode<T>* L, ADnode<T>* R): ADnode<T>("ADPOW",L,R){
  T cv = pow(this->ptrL->curval,this->ptrR->curval);
  this->setValue(cv);
}
template<class T>
ADpow<T>::ADpow(ADnode<T> L,ADnode<T> R): ADnode<T>("ADPOW",&L,&R){
  T cv = pow(this->ptrL->curval,this->ptrR->curval);
  this->setValue(cv);
}

template<class T>
T ADpow<T>::fn(vector<T> x){
  return pow(this->ptrL->fn(x), this->ptrR->fn(x));
}

template<class T>
vector<T> ADpow<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  vector<T> dfnR = this->ptrR->dfn(x);
  T fnL = this->ptrL->fn(x);
  T fnR = this->ptrR->fn(x);
  T pfn = pow(fnL,fnR);
  for(int i = 0; i < x.size(); ++i){
    res[i] = pfn * (dfnR[i] * log(fnL) + fnR * dfnL[i] / fnL);
  }
  return res;      
}


template<class T>
AD<T> pow(const AD<T>& x, const AD<T>& y){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  ADnode<T>* ory = y.getRoot();
  newAD.setRoot(new ADpow<T>(orx,ory));
  return newAD;
}

template<class T>
AD<T> pow(const AD<T>& x, double y){
  AD<T> RAD(y);
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADpow<T>(orx,RAD.getRoot()));
  return newAD;
}


template<class T>
void ADpow<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T fnR = this->ptrR->getValue();
  T wL = w * pow(fnL,fnR - T(1.0)) * fnR;
  T wR = w * pow(fnL,fnR) * log(fnL);
  this->ptrL->bdfn(wL,theta);
  this->ptrR->bdfn(wR,theta);
  return;
}


template<class T>
void ADpow<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  this->ptrR->update(theta);
  T cv = pow(this->ptrL->curval,this->ptrR->curval);
  this->setValue(cv);
  return;
}

