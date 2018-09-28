// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>


// Base class
// template<class T>
// struct AD;

// Base AD graph classes

template<class T>
class ADnode;
template<class T>
struct ADparlist;

// Operations

template<class T>
struct ADparameter;
template<class T>
struct ADconstant;
template<class T>
struct ADsum;
template<class T>
struct ADprod;
template<class T>
struct ADsubtract;

template<class T>
class AD;
