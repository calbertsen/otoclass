namespace my_atomic {


template<class T> int R_finite(T x) { return std::isfinite(asDouble(x)); }
template<class T> int isnan(T x) { return std::isnan(asDouble(x)); }

#undef ML_ERROR
#undef MATHLIB_ERROR
#undef MATHLIB_WARNING
#undef MATHLIB_WARNING2
#undef MATHLIB_WARNING3
#undef MATHLIB_WARNING4
#undef MATHLIB_WARNING5
#undef ML_POSINF
#undef ML_NEGINF
#undef ML_NAN
#undef M_SQRT_2dPI
#undef ISNAN
# define ML_ERROR(x, s) /* nothing */
# define MATHLIB_ERROR(fmt,x) /* nothing */
# define MATHLIB_WARNING(fmt,x) /* nothing */
# define MATHLIB_WARNING2(fmt,x,x2) /* nothing */
# define MATHLIB_WARNING3(fmt,x,x2,x3) /* nothing */
# define MATHLIB_WARNING4(fmt,x,x2,x3,x4) /* nothing */
# define MATHLIB_WARNING5(fmt,x,x2,x3,x4,x5) /* nothing */
#define ML_POSINF	R_PosInf
#define ML_NEGINF	R_NegInf
#define ML_NAN		R_NaN
#define M_SQRT_2dPI	0.797884560802865355879892119869	/* sqrt(2/pi) */
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) */
#define ISNAN(x) (isnan(x)!=0)

#define ML_ERR_return_NAN return R_NaN

#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0) /* 1 */
# define attribute_hidden __attribute__ ((visibility ("hidden")))


  
#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

#define R_D_Cval(p) (lower_tail ? (0.5 - (p) + 0.5) : (p)) /*  1 - p */
  
// Implemented functions //

// lbeta.hpp
template<class Float>
attribute_hidden
Float lbeta_raw(Float a, Float b);

// pt.hpp
template<class Float>
attribute_hidden
Float pt_raw(Float x, Float n, Float lower_tail_, Float log_p_);

// pnorm.hpp
template<class Float>
attribute_hidden
void pnorm_both_raw(Float x, Float *cum, Float *ccum, int i_tail, int log_p);

template<class Float>
attribute_hidden
Float pnorm5_raw(Float x, Float mu, Float sigma, Float lower_tail_, Float log_p_);

}
