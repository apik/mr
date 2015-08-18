/* Derivative of U function with respect to s. */

#include "internal.h"

/* **************************************************************** */

TSIL_COMPLEX TSIL_dUds (TSIL_UTYPE U, TSIL_COMPLEX s)
{
  TSIL_COMPLEX temp, num_th, num_ps, result;
  
  temp = *(U.sval) + *(U.tval[1]) * U.cT2 + *(U.tval[2]) * U.cT3
         + U.con - 0.125L*s;

  num_th = (U.value)*U.cU_th + *(U.tval[0])*U.cT1_th + temp*U.cS_th;
  num_ps = (U.value)*U.cU_ps + *(U.tval[0])*U.cT1_ps + temp*U.cS_ps;

  result = num_th/(s - U.den_th) + num_ps/(s - U.den_ps);
  result /= s;

  return result;
}
