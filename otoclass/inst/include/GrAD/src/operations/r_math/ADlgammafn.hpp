// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>



template<class T>
struct ADlgammafn : public ADnode<T>{


  ADlgammafn(ADnode<T>* L);
  ADlgammafn(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADlgammafn<T>::ADlgammafn(ADnode<T>* L): ADnode<T>("ADLGAMMAFN",L){
  T cv = lgammafn(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADlgammafn<T>::ADlgammafn(ADnode<T> L): ADnode<T>("ADLGAMMAFN",&L){
  T cv = lgammafn(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADlgammafn<T>::fn(vector<T> x){
  return lgammafn(this->ptrL->fn(x));
}

template<class T>
vector<T> ADlgammafn<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * psigamma(fnL,0.0);
  }
  return res;      
}


template<class T>
AD<T> lgammafn(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADlgammafn<T>(orx));
  return newAD;
}



template<class T>
void ADlgammafn<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * psigamma(fnL,0.0);
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADlgammafn<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = lgammafn(this->ptrL->curval);
  this->setValue(cv);
  return;
}

