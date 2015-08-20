/* Derivative of B function with respect to s used for Runge-Kutta.
   This one relies on the TSIL_DATA struct. See the function TSIL_dBds
   in analyticAB.c for a stand-alone version. */

#include "internal.h"

/* **************************************************************** */

TSIL_COMPLEX TSIL_dBds_rk (TSIL_BTYPE B, TSIL_COMPLEX s)
{
  return ((B.B_cB[0] * (B.value) - 0.5L*s + B.B_c[0])/(s - B.B_den[0]) 
        + (B.B_cB[1] * (B.value) - 0.5L*s + B.B_c[1])/(s - B.B_den[1]))/s;
}
