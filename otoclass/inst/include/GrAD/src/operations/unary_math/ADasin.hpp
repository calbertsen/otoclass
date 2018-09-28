// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADasin : public ADnode<T>{


  ADasin(ADnode<T>* L);
  ADasin(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};

template<class T>
ADasin<T>::ADasin(ADnode<T>* L): ADnode<T>("ADASIN",L){
  T cv = asin(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADasin<T>::ADasin(ADnode<T> L): ADnode<T>("ADASIN",&L){
  T cv = asin(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADasin<T>::fn(vector<T> x){
  return asin(this->ptrL->fn(x));
}

template<class T>
vector<T> ADasin<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] / sqrt(-fnL*fnL+T(1.0));
  }
  return res;      
}

template<class T>
AD<T> asin(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADasin<T>(orx));
  return newAD;
}


template<class T>
void ADasin<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w / sqrt(-fnL*fnL + T(1.0));
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADasin<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = asin(this->ptrL->curval);
  this->setValue(cv);
  return;
}

