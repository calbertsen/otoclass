// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADacosh : public ADnode<T>{


  ADacosh(ADnode<T>* L);
  ADacosh(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);


};

template<class T>
ADacosh<T>::ADacosh(ADnode<T>* L): ADnode<T>("ADACOSH",L){
  T cv = acosh(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADacosh<T>::ADacosh(ADnode<T> L): ADnode<T>("ADACOSH",&L){
  T cv = acosh(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADacosh<T>::fn(vector<T> x){
  return acosh(this->ptrL->fn(x));
}

template<class T>
vector<T> ADacosh<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] / (sqrt(fnL+T(1.0)) * sqrt(fnL-T(1.0)));
  }
  return res;      
}

template<class T>
AD<T> acosh(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADacosh<T>(orx));
  return newAD;
}


template<class T>
void ADacosh<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w / (sqrt(fnL + T(1.0)) * sqrt(fnL - T(1.0)));
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADacosh<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = acosh(this->ptrL->curval);
  this->setValue(cv);
  return;
}

