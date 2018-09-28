// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADerf : public ADnode<T>{


  ADerf(ADnode<T>* L);
  ADerf(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);


};


template<class T>
ADerf<T>::ADerf(ADnode<T>* L): ADnode<T>("ADERF",L){
  T cv = erf(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADerf<T>::ADerf(ADnode<T> L): ADnode<T>("ADERF",&L){
  T cv = erf(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADerf<T>::fn(vector<T> x){
  return erf(this->ptrL->fn(x));
}

template<class T>
vector<T> ADerf<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * exp(-fnL*fnL) * T(2.0/sqrt(M_PI));
  }
  return res;      
}


template<class T>
AD<T> erf(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADerf<T>(orx));
  return newAD;
}


template<class T>
void ADerf<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * exp(-fnL*fnL) * T(2.0/sqrt(M_PI));
  this->ptrL->bdfn(wL,theta);
  return;
}

template<class T>
void ADerf<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = erf(this->ptrL->curval);
  this->setValue(cv);
  return;
}

