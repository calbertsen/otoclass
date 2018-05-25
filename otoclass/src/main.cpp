#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS

#define R_NO_REMAP
#include <Rinternals.h>

#include "../inst/include/convert.hpp"
#include "../inst/include/convol.hpp"
#include "../inst/include/efd.hpp"
#include "../inst/include/knn.hpp"

#define install Rf_install
#define findVar Rf_findVar

#define WITH_LIBTMB
#include <TMB.hpp>

extern "C" {

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
    
    CALLDEF(convol2d,2),
    CALLDEF(knn,5),
    CALLDEF(efd2coordSEXP,2),
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