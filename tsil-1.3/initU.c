/* Routines for setup and initial values of U-type functions */

#include "internal.h"

/* **************************************************************** */

void TSIL_ConstructU (TSIL_UTYPE *u,
		 int n,
		 TSIL_REAL z,
		 TSIL_REAL x,
		 TSIL_REAL y,
		 TSIL_REAL v,
		 TSIL_REAL qq)
{
  u->which  = n;
  u->arg[0] = z;
  u->arg[1] = x;
  u->arg[2] = y;
  u->arg[3] = v;

  /*=== Precompute evolution coefficients: ===*/

  /* If second argument vanishes, we will run U=0 exactly, and fix up
     the result afterwards with CorrectUs. */
  if (x < TSIL_TOL) 
  {
    /* Denominator poles (safely avoids s = pole because of rescaling) */
    u->den_th = 2.0L;
    u->den_ps = 2.0L;

    /* Numerator factors: */
    u->cU_th = 0.0L;
    u->cU_ps = 0.0L;
    u->cS_th = 0.0L;
    u->cS_ps = 0.0L;
    u->cT1_th = 0.0L;
    u->cT1_ps = 0.0L;
    u->cT2 = 0.0L;
    u->cT3 = 0.0L;
    u->con = 0.0L;

    return;
  }

  /* Denominator factors */
  u->den_th = TSIL_Th2(z, x);
  u->den_ps = TSIL_Ps2(z, x);

  /* Numerator factors: */
  u->cU_th = u->den_th/2.0L;
  u->cU_ps = u->den_ps/2.0L;

  u->cS_th = -TSIL_SQRT(z/x) - 1.0L;
  u->cS_ps = TSIL_SQRT(z/x) - 1.0L;

  u->cT1_th = (u->cS_th) * (z + 0.5L*TSIL_SQRT(z*x));
  u->cT1_ps = (u->cS_ps) * (z - 0.5L*TSIL_SQRT(z*x));

  u->cT2 = v/2.0L;
  u->cT3 = y/2.0L;

  u->con = 0.5L*(z - TSIL_A(z, qq) + y - TSIL_A(y, qq) + v - TSIL_A(v, qq) -TSIL_I2(x, y, v, qq));

  return;
}
