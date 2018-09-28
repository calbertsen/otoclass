// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADfabs : public ADnode<T>{


  ADfabs(ADnode<T>* L);
  ADfabs(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADfabs<T>::ADfabs(ADnode<T>* L): ADnode<T>("ADFABS",L){
  T cv = fabs(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADfabs<T>::ADfabs(ADnode<T> L): ADnode<T>("ADFABS",&L){
  T cv = fabs(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADfabs<T>::fn(vector<T> x){
  return fabs(this->ptrL->fn(x));
}

template<class T>
vector<T> ADfabs<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i]*sign(fnL);
  }
  return res;      
}

template<class T>
AD<T> fabs(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADfabs<T>(orx));
  return newAD;
}


template<class T>
void ADfabs<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * sign(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}

template<class T>
void ADfabs<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = fabs(this->ptrL->curval);
  this->setValue(cv);
  return;
}

