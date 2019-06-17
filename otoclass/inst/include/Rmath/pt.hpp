/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 2000-2007   The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/*
 * Modified 2019 Christoffer Moesgaard Albertsen
 */

namespace my_atomic {
  
  template<class Float>
  attribute_hidden
  Float pt_raw(Float x, Float n, Float lower_tail_, Float log_p_)
  {
    /* return  P[ T <= x ]	where
     * T ~ t_{n}  (t distrib. with n degrees of freedom).
     *	--> ./pnt.c for NON-central
     */
     int lower_tail = (int)trunc(lower_tail_);
int log_p = (int)trunc(log_p_);
    Float val, nx;
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(n))
      return x + n;
#endif
    if (n <= 0.0) ML_ERR_return_NAN;

    if(!R_FINITE(x))
      return (x < 0) ? R_DT_0 : R_DT_1;
    if(!R_FINITE(n))
      return pnorm5_raw((Float)x, (Float)0.0, (Float)1.0, lower_tail_, log_p_);

#ifdef R_version_le_260
    if (n > 4e5) { /*-- Fixme(?): test should depend on `n' AND `x' ! */
      /* Approx. from	 Abramowitz & Stegun 26.7.8 (p.949) */
      val = 1./(4.*n);
      return pnorm5_raw((Float)(x*(1. - val)/sqrt(1. + x*x*2.*val)), (Float)(0.0), (Float)(1.0),
		   lower_tail_, log_p_);
    }
#endif

    nx = 1 + (x/n)*x;
    /* FIXME: This test is probably losing rather than gaining precision,
     * now that pbeta(*, log_p = TRUE) is much better.
     * Note however that a version of this test *is* needed for x*x > D_MAX */
    if(nx > 1e100) { /* <==>  x*x > 1e100 * n  */
      /* Danger of underflow. So use Abramowitz & Stegun 26.5.4
	 pbeta(z, a, b) ~ z^a(1-z)^b / aB(a,b) ~ z^a / aB(a,b),
	 with z = 1/nx,  a = n/2,  b= 1/2 :
      */
      Float lval;
      lval = -0.5*n*(2*log(fabs(x)) - log(n))
	- lbeta_raw((Float)(0.5*n), (Float)0.5) - log(0.5*n);
      val = log_p ? lval : exp(lval);
    } else {
      val = (n > x * x)
	? atomic::toms708::pbeta((Float)(x * x / (n + x * x)), (Float)(0.5), (Float)(n / 2.), /*lower_tail*/0, log_p)
	: atomic::toms708::pbeta((Float)(1. / nx), (Float)(n / 2.), (Float)(0.5), /*lower_tail*/1, log_p);
    }

    /* Use "1 - v"  if	lower_tail  and	 x > 0 (but not both):*/
    if(x <= 0.)
      lower_tail = !lower_tail;

    if(log_p) {
      if(lower_tail) return log1p(-0.5*exp(val));
      else return val - M_LN2; /* = log(.5* pbeta(....)) */
    }
    else {
      val /= 2.;
      return R_D_Cval(val);
    }
  }

  template<class Float>
  Float pt1(Float x, Float n, Float lower_tail, Float log_p)
  {
    return pt_raw(x,n,lower_tail,log_p);
  }

  TMB_BIND_ATOMIC(pt2,1100,pt1(x[0], x[1], x[2], x[3]))

}
