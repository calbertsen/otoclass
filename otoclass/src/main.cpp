#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#define R_NO_REMAP
#include <Rinternals.h>

#include "../inst/include/biascorrect.hpp"
#include "../inst/include/convert.hpp"
#include "../inst/include/convol.hpp"
#include "../inst/include/efd.hpp"
#include "../inst/include/knn.hpp"
#include "../inst/include/polygon.hpp"
#include "../inst/include/floodfill.hpp"

#define mkString Rf_mkString
#define mkChar Rf_mkChar
#define allocVector Rf_allocVector
#define ScalarInteger Rf_ScalarInteger
#define ScalarLogical Rf_ScalarLogical
#define isNull Rf_isNull
#define install Rf_install
#define findVar Rf_findVar

#define WITH_LIBTMB
#include <TMB.hpp>

extern "C" {


  SEXP bcsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcdsplineR(SEXP x, SEXP knots, SEXP pars);
  SEXP ibcisplineR(SEXP x, SEXP knots, SEXP pars);
  
#include <R_ext/Rdynload.h>
#define CALLDEF(name,n) {#name, (DL_FUNC) &name, n}
 
  static const
  R_CallMethodDef callMethods[] = {
#ifdef TMB_CALLDEFS
    TMB_CALLDEFS,
#else
    CALLDEF(FreeADFunObject, 1)
    CALLDEF(MakeADFunObject, 4),
    CALLDEF(InfoADFunObject, 1),
    CALLDEF(EvalADFunObject, 3),
    CALLDEF(MakeDoubleFunObject, 3),
    CALLDEF(EvalDoubleFunObject, 3),
    CALLDEF(getParameterOrder, 3),
    CALLDEF(MakeADGradObject, 3),
    CALLDEF(MakeADHessObject2, 4),
    CALLDEF(usingAtomics, 0),
    CALLDEF(TMBconfig, 2),
#endif

    CALLDEF(biascorrect_logistic,1),
    CALLDEF(biascorrect_logistic_gradient,1),
    CALLDEF(biascorrect_objective,3),
    CALLDEF(biascorrect_gradient,3),
    CALLDEF(convol2d,2),
    CALLDEF(knn,5),
    CALLDEF(polygon_area,2),
    CALLDEF(point_in_polygon,4),
    CALLDEF(efd2coordSEXP,4),
    CALLDEF(efd,3),
    CALLDEF(scanlineFill,5),
    CALLDEF(bcsplineR,3),
    CALLDEF(ibcsplineR,3),
    CALLDEF(ibcdsplineR,3),
    CALLDEF(ibcisplineR,3),
    {NULL,NULL,0}
  };

  void R_init_otoclass(DllInfo *info)
  {
    /* Register the .C and .Call routines.
       No .Fortran() or .External() routines,
       so pass those arrays as NULL.
    */
    R_registerRoutines(info,
  		       NULL, callMethods,
  		       NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info,FALSE);
  }

  void R_unload_otoclass(DllInfo *info)
  {

  }
}
