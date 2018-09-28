// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADtan : public ADnode<T>{


  ADtan(ADnode<T>* L);
  ADtan(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};

template<class T>
ADtan<T>::ADtan(ADnode<T>* L): ADnode<T>("ADTAN",L){
  T cv = tan(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADtan<T>::ADtan(ADnode<T> L): ADnode<T>("ADTAN",&L){
  T cv = tan(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADtan<T>::fn(vector<T> x){
  return tan(this->ptrL->fn(x));
}

template<class T>
vector<T> ADtan<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * (T(1.0) + tan(fnL) * tan(fnL));
  }
  return res;      
}


template<class T>
AD<T> tan(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADtan<T>(orx));
  return newAD;
}


template<class T>
void ADtan<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * (T(1.0) + tan(fnL) * tan(fnL));
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADtan<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = tan(this->ptrL->curval);
  this->setValue(cv);
  return;
}

