// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>




double sign(double x){
  if(x < 0.0){
    return -1.0;
  }else if(x > 0.0){
    return 1.0;
  }else{
    return 0.0;
  }
}

float sign(float x){
  if(x < 0.0){
    return -1.0;
  }else if(x > 0.0){
    return 1.0;
  }else{
    return 0.0;
  }
}

long double sign(long double x){
  if(x < 0.0){
    return -1.0;
  }else if(x > 0.0){
    return 1.0;
  }else{
    return 0.0;
  }
}


template<class T>
struct ADsign : public ADnode<T>{


  ADsign(ADnode<T>* L);
  ADsign(ADnode<T> L);

  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADsign<T>::ADsign(ADnode<T>* L): ADnode<T>("ADSIGN",L){
  T cv = sign(this->ptrL->curval);
  this->setValue(cv);
}
template<class T>
ADsign<T>::ADsign(ADnode<T> L): ADnode<T>("ADSIGN",&L){
  T cv = sign(this->ptrL->curval);
  this->setValue(cv);
}

template<class T>
T ADsign<T>::fn(vector<T> x){
  return sign(this->ptrL->fn(x));
}

template<class T>
vector<T> ADsign<T>::dfn(vector<T> x){
  vector<T> res(x.size());
  vector<T> dfnL = this->ptrL->dfn(x);
  T fnL = this->ptrL->fn(x);
  for(int i = 0; i < x.size(); ++i){
    res[i] = T(0.0);
  }
  return res;      
}

template<class T>
AD<T> sign(const AD<T>& x){
  AD<T> newAD;
  ADnode<T>* orx = x.getRoot();
  newAD.setRoot(new ADsign<T>(orx));
  return newAD;
}


template<class T>
void ADsign<T>::bdfn(T w, vector<T>& theta){
  // Derivative is zero; no need to continue.
  // T fnL = this->ptrL->getValue();
  // T wL = 0.0;
  // this->ptrL->bdfn(wL,theta);
  return;
}


template<class T>
void ADsign<T>::update(vector<T>& theta){
  this->ptrL->update(theta);
  T cv = sign(this->ptrL->curval);
  this->setValue(cv);
  return;
}
