// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADcos : public ADnode<T>{


  ADcos(ADnode<T>* L);
  ADcos(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};

template<class T>
ADcos<T>::ADcos(ADnode<T>* L): ADnode<T>("ADCOS",L){
  T cv = cos(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADcos<T>::ADcos(ADnode<T> L): ADnode<T>("ADCOS",&L){
  T cv = cos(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADcos<T>::fn(vector<T> x){
  return cos(this->ptrL->fn(x));
}

template<class T>
vector<T> ADcos<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = - dfnL[i] * sin(fnL);
  }
  return res;      
}


template<class T>
AD<T> cos(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADcos<T>(orx));
  return newAD;
}


template<class T>
void ADcos<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = -w * sin(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADcos<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = cos(this->ptrL->curval);
  this->setValue(cv);
  return;
}
