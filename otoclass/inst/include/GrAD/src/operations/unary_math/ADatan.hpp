// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADatan : public ADnode<T>{


  ADatan(ADnode<T>* L);
  ADatan(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);


};

template<class T>
ADatan<T>::ADatan(ADnode<T>* L): ADnode<T>("ADATAN",L){
  T cv = atan(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADatan<T>::ADatan(ADnode<T> L): ADnode<T>("ADATAN",&L){
  T cv = atan(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADatan<T>::fn(vector<T> x){
  return atan(this->ptrL->fn(x));
}

template<class T>
vector<T> ADatan<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] / (fnL*fnL+T(1.0));
  }
  return res;      
}


template<class T>
AD<T> atan(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADatan<T>(orx));
  return newAD;
}


template<class T>
void ADatan<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w / (fnL*fnL+T(1.0));
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADatan<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = atan(this->ptrL->curval);
  this->setValue(cv);
  return;
}

