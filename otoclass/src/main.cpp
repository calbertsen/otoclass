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

#define install Rf_install
#define findVar Rf_findVar

#define WITH_LIBTMB
#include <TMB.hpp>

extern "C" {

  
SEXP logspace_add_p (SEXP logx, SEXP logy, SEXP p);
SEXP logdkwmix(SEXP x, SEXP p, SEXP df, SEXP a, SEXP b);
SEXP logpkwmix(SEXP x, SEXP p, SEXP df, SEXP a, SEXP b);
SEXP logdkwmixE0(SEXP x, SEXP p, SEXP df, SEXP a, SEXP b);
SEXP logpkwmixE0(SEXP x, SEXP p, SEXP df, SEXP a, SEXP b);
  

#include <R_ext/Rdynload.h>
#define CALLDEF(name,n) {#name, (DL_FUNC) &name, n}
 
  static const
  R_CallMethodDef callMethods[] = {
#ifdef TMB_CALLDEFS
    TMB_CALLDEFS,
#else
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
    CALLDEF(polygon_centroid,2),
    CALLDEF(efd2coordSEXP,4),
    CALLDEF(efd,3),
    CALLDEF(scanlineFill,5),

    CALLDEF(logspace_add_p,3),
    CALLDEF(logdkwmix,5),
    CALLDEF(logpkwmix,5),
    CALLDEF(logdkwmixE0,5),
    CALLDEF(logpkwmixE0,5),
    
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
