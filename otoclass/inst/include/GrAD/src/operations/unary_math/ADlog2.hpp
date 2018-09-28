// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADlog2 : public ADnode<T>{


  ADlog2(ADnode<T>* L);
  ADlog2(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADlog2<T>::ADlog2(ADnode<T>* L): ADnode<T>("ADLOG2",L){
  T cv = log2(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADlog2<T>::ADlog2(ADnode<T> L): ADnode<T>("ADLOG2",&L){
  T cv = log2(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADlog2<T>::fn(vector<T> x){
  return log2(this->ptrL->fn(x));
}

template<class T>
vector<T> ADlog2<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] / (fnL * log(T(2.0)));
  }
  return res;      
}


template<class T>
AD<T> log2(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADlog2<T>(orx));
  return newAD;
}


template<class T>
void ADlog2<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w / (fnL * log(T(2.0)));
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADlog2<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = log2(this->ptrL->curval);
  this->setValue(cv);
  return;
}

