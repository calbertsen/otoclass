// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADabs : public ADnode<T>{


  ADabs(ADnode<T>* L);
  ADabs(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);


};


template<class T>
ADabs<T>::ADabs(ADnode<T>* L): ADnode<T>("ADABS",L){
  T cv = abs(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADabs<T>::ADabs(ADnode<T> L): ADnode<T>("ADABS",&L){
  T cv = abs(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADabs<T>::fn(vector<T> x){
  return abs(this->ptrL->fn(x));
}

template<class T>
vector<T> ADabs<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i]*sign(fnL);
  }
  return res;      
}

template<class T>
AD<T> abs(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADabs<T>(orx));
  return newAD;
}


template<class T>
void ADabs<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * sign(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADabs<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = abs(this->ptrL->curval);
  this->setValue(cv);
  return;
}

