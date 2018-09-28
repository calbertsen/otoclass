// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADatan2 : public ADnode<T>{


  ADatan2(ADnode<T>* L, ADnode<T>* R);
  ADatan2(ADnode<T> L,ADnode<T> R);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};

template<class T>
ADatan2<T>::ADatan2(ADnode<T>* L, ADnode<T>* R): ADnode<T>("ADATAN2",L,R){
  T cv = atan2(this->ptrL->curval,this->ptrR->curval);
  this->setValue(cv);
}
template<class T>
ADatan2<T>::ADatan2(ADnode<T> L,ADnode<T> R): ADnode<T>("ADATAN2",&L,&R){
  T cv = atan2(this->ptrL->curval,this->ptrR->curval);
  this->setValue(cv);
}

template<class T>
T ADatan2<T>::fn(vector<T> x){
  return atan2(this->ptrL->fn(x), this->ptrR->fn(x));
}

template<class T>
vector<T> ADatan2<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  vector<T> dfnR = this->ptrR->dfn(x);
  T fnL = this->ptrL->fn(x);
  T fnR = this->ptrR->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = (dfnL[i]*fnR - dfnR[i]*fnL) / (fnL*fnL + fnR*fnR);
  }
  return res;      
}


template<class T>
AD<T> atan2(const AD<T>& x, const AD<T>& y){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  ADnode<T>* ory = y.getRoot();
  newAD.setRoot(new ADatan2<T>(orx,ory));
  return newAD;
}


template<class T>
void ADatan2<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T fnR = this->ptrR->getValue();
  T wL = w / (fnR + fnR * (fnL*fnL)/(fnR*fnR)) ;
  T wR = w * (-fnL / (fnR*fnR + fnR*fnR * (fnL*fnL)/(fnR*fnR)));
  this->ptrL->bdfn(wL,theta);
  this->ptrR->bdfn(wR,theta);
  return;
}

template<class T>
void ADatan2<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  this->ptrR->update(theta);
  T cv = atan2(this->ptrL->curval,this->ptrR->curval);
  this->setValue(cv);
  return;
}

