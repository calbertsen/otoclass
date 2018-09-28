// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>



template<class T>
class ADnode {
  friend class AD<T>;

  
private:
  int parents;
  int owners;

public:
  
  ADparlist<T>* grph;

  bool isBaseParam;
  string baseName;
  int whichParam;
  ADnode<T>* ptrL;
  ADnode<T>* ptrR;
  T curval;

  
  virtual T fn(vector<T> x);
  virtual vector<T> dfn(vector<T> x);

  virtual void bdfn(T w, vector<T>& theta);
  virtual void update(vector<T>& theta);

  virtual string JSON();
  
  void setValue(T x0);
  T getValue();

  
protected:
  
  void prune(AD<T>* requester);
  void addOwner(const AD<T>* o);
  void removeOwner(const AD<T>* o);
  
  void addParent(ADnode<T>* p);
  void removeParent(ADnode<T>* p);
  
protected:
  ADnode();
  ADnode(const string& name, int paramNum, ADparlist<T>* graph);
  ADnode(const string& name);
  ADnode(const string& name, ADnode<T>* L, ADnode<T>* R);
  ADnode(const string& name, ADnode<T>* L);
  ADnode(const string& name, ADnode<T>& L, ADnode<T>& R);
  virtual ~ADnode();

};

template<class T>
ADnode<T>::~ADnode(){}

template<class T>
ADnode<T>::ADnode() : parents(0), owners(0), isBaseParam(false), baseName("UNINITIALIZED")
{
  grph = NULL;
  whichParam = -100;
  ptrL = NULL;
  ptrR = NULL;
}

template<class T>
ADnode<T>::ADnode(const string& name, int paramNum,ADparlist<T>* graph) : parents(0), owners(0), isBaseParam(true), baseName(name), whichParam(paramNum)
{ 
  grph = graph;
  if(whichParam == -1)
    throw 101;
  ptrL = NULL;
  ptrR = NULL;
}

template<class T>
ADnode<T>::ADnode(const string& name) : parents(0), owners(0), isBaseParam(false), baseName(name)
{ 
  grph = NULL;
  whichParam = -1;
  ptrL = NULL;
  ptrR = NULL;
}


template<class T>
ADnode<T>::ADnode(const string& name, ADnode<T>* L, ADnode<T>* R) : parents(0), owners(0), isBaseParam(false), baseName(name) {
  L->addParent(this);
  R->addParent(this);
  if((L->grph != R->grph) && L->grph != NULL && R->grph != NULL)
    throw 99; // Graphs must match

  if(L->grph != NULL){
    grph = L->grph;
  }else if(R->grph != NULL){
    grph = R->grph;
  }else{
    grph = NULL;
  }
  ptrL = L;
  ptrR = R;
  whichParam = -100;
}

template<class T>
ADnode<T>::ADnode(const string& name, ADnode<T>* L) : parents(0), owners(0), isBaseParam(false), baseName(name) {
  L->addParent(this);
  grph = L->grph;
  ptrL = L;
  ptrR = NULL;
  whichParam = -100;
}

template<class T>
ADnode<T>::ADnode(const string& name, ADnode<T>& L, ADnode<T>& R) : parents(0), owners(0), isBaseParam(false), baseName(name){
  L.addParent(this);
  R.addParent(this);
  if(L.grph != R.grph)
    throw 99; // Graphs must match

  if(L.grph != NULL){
    grph = L.grph;
  }else if(R.grph != NULL){
    grph = R.grph;
  }else{
    grph = NULL;
  }
  ptrL = &L;
  ptrR = &R;
  whichParam = -100;
}


template<class T>
T ADnode<T>::fn(vector<T> x){
  throw 100;
}

template<class T>
vector<T> ADnode<T>::dfn(vector<T> x){
  throw 100;
}

template<class T>
void ADnode<T>::bdfn(T w, vector<T>& theta){
  throw 100;
}


template<class T>
string ADnode<T>::JSON(){
  if(ptrL != NULL && ptrR != NULL){
    string sl = ptrL->JSON();
    string sr = ptrR->JSON();
    return "{\"name\":\""+baseName+"\", \"left\":"+sl+", \"right\":"+sr+"}";
  }else if(ptrL != NULL && ptrR == NULL){
    string sl = ptrL->JSON();
    return "{\"name\":\""+baseName+"\", \"left\":"+sl+", \"right\":null}";
  }else if(ptrL == NULL && ptrR != NULL){
    string sr = ptrR->JSON();
    return "{\"name\":\""+baseName+"\", \"left\":null, \"right\":"+sr+"}";
  }else{
    return "{\"name\":\""+baseName+"\", \"left\":null, \"right\":null}";
  }

}

template<class T>
void ADnode<T>::update(vector<T>& theta){
  throw 999;
}


template<class T>
void ADnode<T>::setValue(T x0){
  curval = x0;
}

template<class T>
T ADnode<T>::getValue(){
  return curval;
}


template<class T>
void ADnode<T>::addParent(ADnode<T>* p){
  parents++;
  return;
}

template<class T>
void ADnode<T>::removeParent(ADnode<T>* p){
  parents--;
  return;
}


template<class T>
void ADnode<T>::addOwner(const AD<T>* o){
  if(o != NULL)
    owners++;
  return;
}

template<class T>
void ADnode<T>::removeOwner(const AD<T>* o){
  if(o != NULL)
    owners--;
  if(owners < 0)
    throw "OUT OF BOUNDS";
  return;
}


template<class T>
void ADnode<T>::prune(AD<T>* requester){
  removeOwner(requester);
  if(parents == 0 && owners == 0){
    if(ptrL != NULL){
      ptrL->removeParent(this);
      ptrL->prune(NULL);
      ptrL = NULL;
    }
    if(ptrR != NULL){
      ptrR->removeParent(this);
      ptrR->prune(NULL);
      ptrR = NULL;
    }
    if(requester != NULL || !isBaseParam)// && !isBaseParam)
      delete this;
  }
}
