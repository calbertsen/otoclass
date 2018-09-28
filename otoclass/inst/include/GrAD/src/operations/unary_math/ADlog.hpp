// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADlog : public ADnode<T>{


  ADlog(ADnode<T>* L);
  ADlog(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADlog<T>::ADlog(ADnode<T>* L): ADnode<T>("ADLOG",L){
  T cv = log(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADlog<T>::ADlog(ADnode<T> L): ADnode<T>("ADLOG",&L){
  T cv = log(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADlog<T>::fn(vector<T> x){
  return log(this->ptrL->fn(x));
}

template<class T>
vector<T> ADlog<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] / fnL;
  }
  return res;      
}


template<class T>
AD<T> log(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADlog<T>(orx));
  return newAD;
}



template<class T>
void ADlog<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w /fnL;
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADlog<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = log(this->ptrL->curval);
  this->setValue(cv);
  return;
}

