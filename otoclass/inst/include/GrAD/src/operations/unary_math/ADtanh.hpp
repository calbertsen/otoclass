// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADtanh : public ADnode<T>{


  ADtanh(ADnode<T>* L);
  ADtanh(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADtanh<T>::ADtanh(ADnode<T>* L): ADnode<T>("ADTANH",L){
  T cv = tanh(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADtanh<T>::ADtanh(ADnode<T> L): ADnode<T>("ADTANH",&L){
  T cv = tanh(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADtanh<T>::fn(vector<T> x){
  return tanh(this->ptrL->fn(x));
}

template<class T>
vector<T> ADtanh<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * (T(1.0) - tanh(fnL) * tanh(fnL));
  }
  return res;      
}


template<class T>
AD<T> tanh(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADtanh<T>(orx));
  return newAD;
}


template<class T>
void ADtanh<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * (T(1.0) - tanh(fnL) * tanh(fnL));
  this->ptrL->bdfn(wL,theta);
  return;
}

template<class T>
void ADtanh<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = tanh(this->ptrL->curval);
  this->setValue(cv);
  return;
}

