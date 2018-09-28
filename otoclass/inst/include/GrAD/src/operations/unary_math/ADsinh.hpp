// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADsinh : public ADnode<T>{


  ADsinh(ADnode<T>* L);
  ADsinh(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADsinh<T>::ADsinh(ADnode<T>* L): ADnode<T>("ADSINH",L){
  T cv = sinh(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADsinh<T>::ADsinh(ADnode<T> L): ADnode<T>("ADSINH",&L){
  T cv = sinh(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADsinh<T>::fn(vector<T> x){
  return sinh(this->ptrL->fn(x));
}

template<class T>
vector<T> ADsinh<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * cosh(fnL);
  }
  return res;      
}


template<class T>
AD<T> sinh(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADsinh<T>(orx));
  return newAD;
}


template<class T>
void ADsinh<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * cosh(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADsinh<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = sinh(this->ptrL->curval);
  this->setValue(cv);
  return;
}

