// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


template<class T>
struct ADparlist {
  
  int nparams;
  vector<AD<T>* > params;

  ADparlist();

  //ADparlist(vector<T*> par);
    
  int which(AD<T>* par);

  int Independent(AD<T>* par);
  int Independent(AD<T>& par);
  int Independent(vector<AD<T>* >& par);
  int Independent(vector<AD<T> >* par);
  int Independent(vector<AD<T> >& par);

};


template<class T>
int ADparlist<T>::Independent(AD<T>* par){
  int parNum = which(par);
  if(parNum < 0){
    params.push_back(par);
    nparams++;
    par->toParameter("PARAMETER_"+itos(nparams-1),nparams-1,this);
    return 1;
  }else{
    return 0;
  }  
}

template<class T>
int ADparlist<T>::Independent(AD<T>& par){
  return Independent(&par);
}

template<class T>
int ADparlist<T>::Independent(vector<AD<T>* >& par){
  int res = 0;
  for(int i = 0; i < par.size(); ++i){
    res += Independent(par[i]);
  }
  return res;
}


template<class T>
int ADparlist<T>::Independent(vector<AD<T> >* par){
  int res = 0;
  for(int i = 0; i < par->size(); ++i){
    res += Independent(&(par->operator[](i)));
  }
  return res;
}

template<class T>
int ADparlist<T>::Independent(vector<AD<T> >& par){
  int res = 0;
  for(int i = 0; i < par.size(); ++i){
    res += Independent(&(par[i]));
  }
  return res;
}



template<class T>
ADparlist<T>::ADparlist(){
  params = vector<AD<T>* >();
  nparams = 0;
}


// template<class T>
// ADparlist<T>::ADparlist(vector<T*> par){
//   params = par;
//   nparams = par.size();
// }

template<class T>
int ADparlist<T>::which(AD<T>* par){
  int res = -1;
  if(nparams == 0)
    return res;
  int i = 0;
  do {
    if(par == params[i])
      res = i;
    ++i;
  } while(res < 0 && i < nparams);
  return res;
}

// void ADparlist::addParam(string par){
//   int i = which(par);
//   if(i == -1){
//     params.push_back(par);
//     nparams++;
//   }
// }


