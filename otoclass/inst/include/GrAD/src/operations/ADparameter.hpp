// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADparameter : public ADnode<T>{

  ADparameter(const string& name, int paramNum, ADparlist<T>* graph);
  ADparameter(T x0, const string& name, int paramNum, ADparlist<T>* graph);

  
  T fn(vector<T> x);
  vector<T> dfn(vector<T> x);
  void bdfn(T w, vector<T>& theta);
  void update(vector<T>& theta);

};


template<class T>
ADparameter<T>::ADparameter(const string& name, int paramNum, ADparlist<T>* graph) :
  ADnode<T>(name, paramNum, graph){
  this->setValue(T(0.0));
}

template<class T>
ADparameter<T>::ADparameter(T x0, const string& name, int paramNum, ADparlist<T>* graph) :
  ADnode<T>(name, paramNum, graph){
  this->setValue(x0);
}

template<class T>
T ADparameter<T>::fn(vector<T> x){
  if(x.size() != this->grph->nparams)
    throw "Wrong argument size to parameter";
  return x[this->whichParam];
}

template<class T>
vector<T> ADparameter<T>::dfn(vector<T> x){
  if(x.size() != this->grph->nparams)
    throw "Wrong argument size to parameter";
  vector<T> res(x.size());
  for(int i = 0; i < x.size(); ++i){
    if(this->whichParam == i){
      res[i] = T(1.0);
    }else{
      res[i] = T(0.0);
    }
  }
  return res;
}

template<class T>
void ADparameter<T>::bdfn(T w, vector<T>& theta){
  theta[this->whichParam] += w;
  return;
}

template<class T>
void ADparameter<T>::update(vector<T>& theta){
  this->setValue(theta[this->whichParam]);
  return;
}
