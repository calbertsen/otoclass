// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// knn
NumericMatrix knn(NumericMatrix train, IntegerVector group, NumericMatrix test, int kn, int disttype);
RcppExport SEXP otoclass_knn(SEXP trainSEXP, SEXP groupSEXP, SEXP testSEXP, SEXP knSEXP, SEXP disttypeSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type train(trainSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type test(testSEXP );
        Rcpp::traits::input_parameter< int >::type kn(knSEXP );
        Rcpp::traits::input_parameter< int >::type disttype(disttypeSEXP );
        NumericMatrix __result = knn(train, group, test, kn, disttype);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
