// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADlog10 : public ADnode<T>{


  ADlog10(ADnode<T>* L);
  ADlog10(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADlog10<T>::ADlog10(ADnode<T>* L): ADnode<T>("ADLOG10",L){
  T cv = log10(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADlog10<T>::ADlog10(ADnode<T> L): ADnode<T>("ADLOG10",&L){
  T cv = log10(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADlog10<T>::fn(vector<T> x){
  return log10(this->ptrL->fn(x));
}

template<class T>
vector<T> ADlog10<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] / (fnL * log(T(10.0)));
  }
  return res;      
}


template<class T>
AD<T> log10(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADlog10<T>(orx));
  return newAD;
}


template<class T>
void ADlog10<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w / (fnL * log(T(10.0)));
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADlog10<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = log10(this->ptrL->curval);
  this->setValue(cv);
  return;
}
