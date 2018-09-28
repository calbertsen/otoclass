// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADlog1p : public ADnode<T>{


  ADlog1p(ADnode<T>* L);
  ADlog1p(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADlog1p<T>::ADlog1p(ADnode<T>* L): ADnode<T>("ADLOG1P",L){
  T cv = log1p(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADlog1p<T>::ADlog1p(ADnode<T> L): ADnode<T>("ADLOG1P",&L){
  T cv = log1p(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADlog1p<T>::fn(vector<T> x){
  return log1p(this->ptrL->fn(x));
}

template<class T>
vector<T> ADlog1p<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] / (fnL + T(1.0));
  }
  return res;      
}


template<class T>
AD<T> log1p(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADlog1p<T>(orx));
  return newAD;
}


template<class T>
void ADlog1p<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w / (fnL + T(1.0));
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADlog1p<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = log1p(this->ptrL->curval);
  this->setValue(cv);
  return;
}

