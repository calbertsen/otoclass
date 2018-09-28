// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADexp2 : public ADnode<T>{


  ADexp2(ADnode<T>* L);
  ADexp2(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADexp2<T>::ADexp2(ADnode<T>* L): ADnode<T>("ADEXP2",L){
  T cv = exp2(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADexp2<T>::ADexp2(ADnode<T> L): ADnode<T>("ADEXP2",&L){
  T cv = exp2(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADexp2<T>::fn(vector<T> x){
  return exp2(this->ptrL->fn(x));
}

template<class T>
vector<T> ADexp2<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * exp2(fnL) * log(T(2.0));
  }
  return res;      
}


template<class T>
AD<T> exp2(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADexp2<T>(orx));
  return newAD;
}


template<class T>
void ADexp2<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * exp2(fnL) * log(T(2.0));
  this->ptrL->bdfn(wL,theta);
  return;
}

template<class T>
void ADexp2<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = exp2(this->ptrL->curval);
  this->setValue(cv);
  return;
}

