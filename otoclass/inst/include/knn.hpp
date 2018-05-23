#include "convert.hpp"

double distance(doubleVector x, doubleVector y,int type);

MatrixXd knn_work(MatrixXd train, intVector group, MatrixXd test, int kn, int disttype, int ngroup);

extern "C" {
  SEXP knn(SEXP train, SEXP group, SEXP test, SEXP kn, SEXP disttype);
}
