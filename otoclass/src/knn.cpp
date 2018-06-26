#include "../inst/include/knn.hpp"

double distance(doubleVector x, doubleVector y,int type){

  double ans = 0.0;
  doubleVector z = x-y;
  
  switch(type){
  case 1: // Euclidian (L_2)
    // ans = sqrt(sum(pow((x-y),2)));
    ans = sqrt(z.pow(2).sum());
  case 2: // Manhattan (L_1)
    // ans = sum(abs(x-y));
    ans = z.abs().sum();
  case 3: // Chebyshev (L_\infty)
    // ans = max(abs(x-y));
    ans = z.abs().maxCoeff();
  }
  return ans;
}

MatrixXd knn_work(MatrixXd train, intVector group, MatrixXd test, int kn, int disttype, int ngroup){
  MatrixXd pred(test.cols(),ngroup);
  pred.setZero();
  for(int i = 0; i < test.cols(); ++i){ // For each test point
    doubleVector dists(kn);
    dists.setZero();
    dists += HUGE_VAL;
    doubleVector indx(kn);
    indx.setZero();
      
    for(int j = 0; j < train.cols(); ++j){ // Run through train points
      double d = distance(train.col(j),test.col(i),disttype); // Find the distance
      if(d < dists(kn-1)){ // If distance is more than the kn furthest away (so far) it can not be one of the kn nearest neighbours. Otherwise insert in list.
	for (int k = 0; k < kn; ++k){ // Run trough list of neighbours 
	  if (d < dists(k)){// If the training point is closer than the one we are looking at, it should be inserted here
	    for (int k1 = kn-1; k1 > k; --k1){ // All neighbours further away should move one place up in the list
	      dists(k1) = dists(k1 - 1);
	      indx(k1) = indx(k1 - 1);
	    }
	    //Finally insert the new point
	    dists(k) = d;
	    indx(k) = j;
	    break;
	  }
	}
      }
    }
    //Run through neighbours
    for(int j = 0; j < indx.size(); ++j){
      // Add a vote to the group it belongs to
      pred(i,group(indx(j))) += 1.0/(double)kn;
    }
  }

  return pred;

}


extern "C" {
  SEXP knn(SEXP train, SEXP group, SEXP test, SEXP kn, SEXP disttype){
    SEXP uniqueGroups;
    PROTECT(uniqueGroups = Rf_eval( Rf_lang2( Rf_install("unique"), group), R_GlobalEnv));
    MatrixXd res = knn_work(asDoubleMatrix(train), asIntVector(group), asDoubleMatrix(test), asInteger(kn), asInteger(disttype), Rf_length(uniqueGroups));
    UNPROTECT(1);
    return asSEXP(res);
  }
}
