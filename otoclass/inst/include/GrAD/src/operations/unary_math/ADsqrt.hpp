// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADsqrt : public ADnode<T>{


  ADsqrt(ADnode<T>* L);
  ADsqrt(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};

template<class T>
ADsqrt<T>::ADsqrt(ADnode<T>* L): ADnode<T>("ADSQRT",L){
  T cv = sqrt(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADsqrt<T>::ADsqrt(ADnode<T> L): ADnode<T>("ADSQRT",&L){
  T cv = sqrt(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADsqrt<T>::fn(vector<T> x){
  return sqrt(this->ptrL->fn(x));
}

template<class T>
vector<T> ADsqrt<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = - T(0.5) * dfnL[i] / sqrt(fnL);
  }
  return res;      
}



template<class T>
AD<T> sqrt(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADsqrt<T>(orx));
  return newAD;
}


template<class T>
void ADsqrt<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = - w * T(0.5) / sqrt(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADsqrt<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = sqrt(this->ptrL->curval);
  this->setValue(cv);
  return;
}

