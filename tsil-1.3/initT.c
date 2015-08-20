/* Routines for setup and initial values of T-type functions */

#include "internal.h"

/* **************************************************************** */

void TSIL_ConstructT (TSIL_TTYPE *t,
		 int n,
		 TSIL_REAL x,
		 TSIL_REAL y,
		 TSIL_REAL z,
		 TSIL_REAL qq)
{
  TSIL_REAL Ay, Az, Ax, Deltaxyz;
  TSIL_REAL sqrtx, sqrty, sqrtz, alphax, alphay, alphaz;
  int i;

  t->which  = n;
  t->arg[0] = x;
  t->arg[1] = y;
  t->arg[2] = z;

  Ay = TSIL_A(y, qq);
  Az = TSIL_A(z, qq);
  Ax = TSIL_A(x, qq);
  Deltaxyz = TSIL_Delta(x,y,z);

  sqrtx = TSIL_SQRT(x);
  sqrty = TSIL_SQRT(y);
  sqrtz = TSIL_SQRT(z);

  alphax = TSIL_Alpha (x, qq);
  alphay = TSIL_Alpha (y, qq);
  alphaz = TSIL_Alpha (z, qq);

  /* If first argument vanishes, disable evolution by setting everything to 0.*/
  if (x < TSIL_TOL) {
    for (i=0; i<4; i++) {
      t->T_den[i] = 0.0L;
      t->cTS_num[i] = 0.0L;
      t->cTT1_num[i] = 0.0L;
      t->cTT2_num[i] = 0.0L;
      t->cTT3_num[i] = 0.0L;
      t->cT_num[i] = 0.0L;
    }
    t->cTs_num = 0.0L;

    return;
  }

  /* Otherwise precompute evolution factors: */

  /* Denominator factors */
  t->T_den[0] = TSIL_POW( sqrtx + sqrty + sqrtz, 2);
  t->T_den[1] = TSIL_POW(-sqrtx + sqrty + sqrtz, 2);
  t->T_den[2] = TSIL_POW( sqrtx - sqrty + sqrtz, 2);
  t->T_den[3] = TSIL_POW( sqrtx + sqrty - sqrtz, 2);

  /* Numerator factors */
  t->cTS_num[0] = (sqrtx + sqrty + sqrtz)/(2.0L*sqrtx);
  t->cTS_num[1] = (sqrtx - sqrty - sqrtz)/(2.0L*sqrtx);
  t->cTS_num[2] = (sqrtx - sqrty + sqrtz)/(2.0L*sqrtx);
  t->cTS_num[3] = (sqrtx + sqrty - sqrtz)/(2.0L*sqrtx);

  /* cTT1 */ 
  t->cTT1_num[0] = 0.25L*(sqrtx + sqrty + sqrtz)*(2.0L*sqrtx + sqrty + sqrtz);
  t->cTT1_num[1] = 0.25L*(sqrtx - sqrty - sqrtz)*(2.0L*sqrtx - sqrty - sqrtz);
  t->cTT1_num[2] = 0.25L*(sqrtx - sqrty + sqrtz)*(2.0L*sqrtx - sqrty + sqrtz);
  t->cTT1_num[3] = 0.25L*(sqrtx + sqrty - sqrtz)*(2.0L*sqrtx + sqrty - sqrtz);

  /* cTT2 */ 
  t->cTT2_num[0] =  (0.25L*sqrty/sqrtx)*(sqrtx + sqrty + sqrtz)*(sqrtx + 2.0L*sqrty + sqrtz);
  t->cTT2_num[1] = -(0.25L*sqrty/sqrtx)*(sqrtx - sqrty - sqrtz)*(sqrtx - 2.0L*sqrty - sqrtz);
  t->cTT2_num[2] = -(0.25L*sqrty/sqrtx)*(sqrtx - sqrty + sqrtz)*(sqrtx - 2.0L*sqrty + sqrtz);
  t->cTT2_num[3] =  (0.25L*sqrty/sqrtx)*(sqrtx + sqrty - sqrtz)*(sqrtx + 2.0L*sqrty - sqrtz);

  /* cTT3 (cTT2 with y <-> z; note switch in denominators though!) */
  t->cTT3_num[0] =  (0.25L*sqrtz/sqrtx)*(sqrtx + sqrtz + sqrty)*(sqrtx + 2.0L*sqrtz + sqrty);
  t->cTT3_num[1] = -(0.25L*sqrtz/sqrtx)*(sqrtx - sqrtz - sqrty)*(sqrtx - 2.0L*sqrtz - sqrty);
  t->cTT3_num[2] =  (0.25L*sqrtz/sqrtx)*(sqrtx + sqrtz - sqrty)*(sqrtx + 2.0L*sqrtz - sqrty);
  t->cTT3_num[3] = -(0.25L*sqrtz/sqrtx)*(sqrtx - sqrtz + sqrty)*(sqrtx - 2.0L*sqrtz + sqrty);

  /* cT */
  t->cT_num[0] = ((sqrtx + sqrty + sqrtz)/(4.0L*sqrtx))*
                 ( 5.75L*x + 0.75L*(y + z)
                 + 2.5L*(sqrtx*sqrty + sqrtx*sqrtz - sqrty*sqrtz)
                 + alphax*(-5.L*sqrtx - 3.L*sqrty - 3.L*sqrtz)
                 + alphay*(sqrtx - sqrty + sqrtz)
                 + alphaz*(sqrtx + sqrty - sqrtz)
                 - alphax*alphay - alphax*alphaz - alphay*alphaz);

  t->cT_num[1] = ((sqrtx - sqrty - sqrtz)/(4.0L*sqrtx))*
                 ( 0.75L*(x + y + z)
                 + 2.5L*(sqrtx*sqrty + sqrtx*sqrtz - sqrty*sqrtz)
                 + alphax*(-sqrtx - sqrty - sqrtz)
                 + alphay*(-sqrtx - sqrty + sqrtz)
                 + alphaz*(-sqrtx + sqrty - sqrtz) 
                 + alphax*alphay + alphax*alphaz - alphay*alphaz);

  t->cT_num[2] = ((sqrtx - sqrty + sqrtz)/(4.0L*sqrtx))*
                 ( 0.75L*(x + y + z)
                 + 2.5L*(sqrtx*sqrty - sqrtx*sqrtz + sqrty*sqrtz)
                 + alphax*(-sqrtx - sqrty + sqrtz)
                 + alphay*(-sqrtx - sqrty - sqrtz)
                 + alphaz*(sqrtx - sqrty - sqrtz)
                 + alphax*alphay - alphax*alphaz + alphay*alphaz);

  t->cT_num[3] = ((sqrtx + sqrty - sqrtz)/(4.0L*sqrtx))*
                 ( 0.75L*(x + y + z) 
                 + 2.5L*(-sqrtx*sqrty + sqrtx*sqrtz + sqrty*sqrtz)
                 + alphax*(-sqrtx + sqrty - sqrtz)
                 + alphay*(sqrtx - sqrty - sqrtz)
                 + alphaz*(-sqrtx - sqrty - sqrtz)
                 - alphax*alphay + alphax*alphaz + alphay*alphaz);

  /* cTs */
  t->cTs_num = alphax/sqrtx - 1.25L;

  return;
}
