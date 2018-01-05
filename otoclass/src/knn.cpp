#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;
  
double distance(NumericVector x, NumericVector y,int type){

  double ans;

  switch(type){
  case 1: // Euclidian (L_2)
    ans = sqrt(sum(pow((x-y),2)));
  case 2: // Manhattan (L_1)
    ans = sum(abs(x-y));
  case 3: // Chebyshev (L_\infty)
    ans = max(abs(x-y));
  }
  return ans;
}

// [[Rcpp::export]]
NumericMatrix knn(NumericMatrix train, IntegerVector group, NumericMatrix test, int kn, int disttype){

  int ngroup = sort_unique(group).size();

  NumericMatrix pred(test.nrow(),ngroup);
  for(int i = 0; i < test.nrow(); ++i){ // For each test point
    NumericVector dists(kn);
    dists = dists + 1000000000.0;
    IntegerVector indx(kn);

      
    for(int j = 0; j < train.nrow(); ++j){ // Run through train points
      double d = distance(train.row(j),test.row(i),disttype); // Find the distance
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
