// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>

template<class T>
class AD {

private:
  
  ADnode<T>* root;

public:
  
  AD();
  AD(const AD<T>& x);
  // AD(T x0, string name, ADparlist<T>* graph);
  AD(T x0);
  ~AD(){
    if(root != NULL)
      root->prune(this);
  }

  void toParameter(const string& name, int paramNum, ADparlist<T>* graph);
  
  T forward_fn(vector<T> x);
  vector<T> forward_gr(vector<T> x);

  T fn();
  vector<T> gr();
  vector<T>& gr(vector<T>& x);
  string JSON();
  void update(vector<T>& x);

  
  void setRoot(ADnode<T>* r){
    //if(root != r){
      ADnode<T>* oldRoot = root;
      root = r;
      if(root != NULL)
     	root->addOwner(this);
      if(oldRoot != NULL)
     	oldRoot->prune(this);
      //}
    return;
  };

  ADnode<T>* getRoot() const{
    return root;
  }

  AD<T>& operator=(const AD<T>& other){
    ADnode<T>* oro = other.getRoot();
    setRoot(oro);
    return *this;
  };
  
  AD<T> operator+(const AD<T>& other) const;
  AD<T>& operator+=(const AD<T>& other);

  AD<T> operator-() const;
  AD<T> operator-(const AD<T>& other) const;
  AD<T>& operator-=(const AD<T>& other);

  AD<T> operator*(const AD<T>& other) const;
  AD<T>& operator*=(const AD<T>& other);

  AD<T> operator/(const AD<T>& other) const;
  AD<T>& operator/=(const AD<T>& other);

  // Friend declarations are not needed while root is public.
  // // Unary math
  // template<class U>
  // friend AD<U> sign(const AD<U>& x);
  // template<class U>
  // friend AD<U> abs(const AD<U>& x);
  // template<class U>
  // friend AD<U> fabs(const AD<U>& x);
  // template<class U>
  // friend AD<U> sqrt(const AD<U>& x);
  // template<class U>
  // friend AD<U> acos(const AD<U>& x);
  // template<class U>
  // friend AD<U> acosh(const AD<U>& x);
  // template<class U>
  // friend AD<U> asin(const AD<U>& x);
  // template<class U>
  // friend AD<U> asinh(const AD<U>& x);
  // template<class U>
  // friend AD<U> atan(const AD<U>& x);
  // template<class U>
  // friend AD<U> atanh(const AD<U>& x);
  // template<class U>
  // friend AD<U> cos(const AD<U>& x);
  // template<class U>
  // friend AD<U> cosh(const AD<U>& x);
  // template<class U>
  // friend AD<U> erf(const AD<U>& x);
  // template<class U>
  // friend AD<U> exp(const AD<U>& x);
  // template<class U>
  // friend AD<U> expm1(const AD<U>& x);
  // template<class U>
  // friend AD<U> log(const AD<U>& x);
  // template<class U>
  // friend AD<U> log1p(const AD<U>& x);
  // template<class U>
  // friend AD<U> log10(const AD<U>& x);
  // template<class U>
  // friend AD<U> sin(const AD<U>& x);
  // template<class U>
  // friend AD<U> sinh(const AD<U>& x);
  // template<class U>
  // friend AD<U> tan(const AD<U>& x);
  // template<class U>
  // friend AD<U> tanh(const AD<U>& x);

  // // Binary math
  // template<class U>
  // friend AD<U> atan2(const AD<U>& x,const AD<U>& y);
  // template<class U>
  // friend AD<U> pow(const AD<U>& x,const AD<U>& y);

};




template<class T>
AD<T>::AD(){
  root = NULL;
}

template<class T>
AD<T>::AD(const AD<T>& x){
  root = NULL;
  ADnode<T>* newRoot = x.getRoot();
  setRoot(newRoot);
}

// template<class T>
// AD<T>::AD(T x0, string ADparlist<T>* graph){
//   int i = graph->which();
//   if(i < 0){
//     root = new ADconstant<T>(x0);
//   }else{
//     root = new ADparameter<T>(x0,name,graph);
//   }
// }

template<class T>
void AD<T>::toParameter(const string& name, int paramNum, ADparlist<T>* graph){
  T x0 = this->root->curval;
  setRoot(new ADparameter<T>(x0,name,paramNum,graph));
  return;
}

template<class T>
AD<T>::AD(T x0){
  root = NULL;
  setRoot(new ADconstant<T>(x0));
  //root->isBaseParam = true; // Quick fix
  return;
}

template<class T>
T AD<T>::forward_fn(vector<T> x){
  return this->root->fn(x);
}

template<class T>
vector<T> AD<T>::forward_gr(vector<T> x){
  return this->root->dfn(x);
}

template<class T>
T AD<T>::fn(){
  return this->root->getValue();
}

template<class T>
vector<T> AD<T>::gr(){
  int np = this->root->grph->nparams;
  vector<T> x(np);
  for(int i = 0; i < np; ++i) x[i] = T(0.0);
  this->root->bdfn(1.0,x);
  return x;
}

template<class T>
vector<T>& AD<T>::gr(vector<T>& x){
  this->root->bdfn(1.0,x);
  return x;
}

template<class T>
string AD<T>::JSON(){
  return this->root->JSON();
}


template<class T>
AD<T> AD<T>::operator+(const AD<T>& other) const{
  AD<T> newAD = AD();
  ADnode<T>* oro = other.getRoot();
  ADnode<T>* ort = getRoot();
  newAD.setRoot(new ADsum<T>(ort,oro));
  return newAD;
}


template<class T>
AD<T>& AD<T>::operator+=(const AD<T>& other){
  ADnode<T>* oro = other.getRoot();
  setRoot(new ADsum<T>(getRoot(),oro));
  return *this;
}

template<class T>
AD<T> AD<T>::operator-() const{
  AD<T> newAD = AD();
  ADnode<T>* ort = getRoot();
  newAD.setRoot(new ADusubtract<T>(ort));
  return newAD;
}


template<class T>
AD<T> AD<T>::operator-(const AD<T>& other) const{
  AD<T> newAD = AD();
  ADnode<T>* oro = other.getRoot();
  ADnode<T>* ort = getRoot();
  newAD.setRoot(new ADsubtract<T>(ort,oro));
  return newAD;
}

template<class T>
AD<T>& AD<T>::operator-=(const AD<T>& other){
  ADnode<T>* oro = other.getRoot();
  setRoot(new ADsubtract<T>(getRoot(),oro));
  return *this;
}

template<class T>
AD<T> AD<T>::operator*(const AD<T>& other) const{
  AD<T> newAD = AD();
  ADnode<T>* oro = other.getRoot();
  ADnode<T>* ort = getRoot();
  newAD.setRoot(new ADprod<T>(ort,oro));
  return newAD;
}

template<class T>
AD<T>& AD<T>::operator*=(const AD<T>& other){
  ADnode<T>* oro = other.getRoot();
  setRoot(new ADprod<T>(getRoot(),oro));
  return *this;
}


template<class T>
AD<T> AD<T>::operator/(const AD<T>& other) const{
  ADnode<T>* oro = other.getRoot();
  AD<T> newAD = AD();
  ADnode<T>* ort = getRoot();
  newAD.setRoot(new ADdiv<T>(ort,oro));
  return newAD;
}

template<class T>
AD<T>& AD<T>::operator/=(const AD<T>& other){
  ADnode<T>* oro = other.getRoot();
  setRoot(new ADdiv<T>(getRoot(),oro));
  return *this;
}

template<class T>
void AD<T>::update(vector<T>& x){
  root->update(x);
}
