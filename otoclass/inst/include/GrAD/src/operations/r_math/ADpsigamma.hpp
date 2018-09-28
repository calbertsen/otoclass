// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>



template<class T>
struct ADpsigamma : public ADnode<T>{

  double p;

  ADpsigamma(ADnode<T>* L, double p_);
  ADpsigamma(ADnode<T> L, double p_);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADpsigamma<T>::ADpsigamma(ADnode<T>* L, double p_): ADnode<T>("ADPSIGAMMA",L), p(p_){
  T cv = psigamma(this->ptrL->curval,p);
  this->setValue(cv);
}
template<class T>
ADpsigamma<T>::ADpsigamma(ADnode<T> L, double p_): ADnode<T>("ADPSIGAMMA",&L), p(p_){
  T cv = psigamma(this->ptrL->curval,p);
  this->setValue(cv);
}

template<class T>
T ADpsigamma<T>::fn(vector<T> x){
  return psigamma(this->ptrL->fn(x),p);
}

template<class T>
vector<T> ADpsigamma<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * psigamma(fnL,p+1.0);
  }
  return res;      
}


template<class T>
AD<T> psigamma(const AD<T>& x, double p){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADpsigamma<T>(orx,p));
  return newAD;
}



template<class T>
void ADpsigamma<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * psigamma(fnL,p+1.0);
  this->ptrL->bdfn(wL,theta);
  return;
}

template<class T>
void ADpsigamma<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = psigamma(this->ptrL->curval,p);
  this->setValue(cv);
  return;
}

