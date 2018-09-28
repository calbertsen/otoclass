// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>



template<class T>
struct ADgammafn : public ADnode<T>{


  ADgammafn(ADnode<T>* L);
  ADgammafn(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADgammafn<T>::ADgammafn(ADnode<T>* L): ADnode<T>("ADGAMMAFN",L){
  T cv = gammafn(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADgammafn<T>::ADgammafn(ADnode<T> L): ADnode<T>("ADGAMMAFN",&L){
  T cv = gammafn(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADgammafn<T>::fn(vector<T> x){
  return gammafn(this->ptrL->fn(x));
}

template<class T>
vector<T> ADgammafn<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = dfnL[i] * psigamma(fnL,0.0)*gammafn(fnL);
  }
  return res;      
}


template<class T>
AD<T> gammafn(AD<T> x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADgammafn<T>(orx));
  return newAD;
}



template<class T>
void ADgammafn<T>::bdfn(T w, vector<T>& theta){
  T fnL = this->ptrL->getValue();
  T wL = w * psigamma(fnL,0.0) * gammafn(fnL);
  this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADgammafn<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = gammafn(this->ptrL->curval);
  this->setValue(cv);
  return;
}

