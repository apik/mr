/* Derivative of S function with respect to s. */

#include "internal.h"

/* ******************************************************************* */

TSIL_COMPLEX TSIL_dSds (TSIL_STYPE S, TSIL_COMPLEX s)
{
  return (S.value + 
	  S.arg[0] * *(S.tval[0]) + 
	  S.arg[1] * *(S.tval[1]) + 
	  S.arg[2] * *(S.tval[2]) + S.S_c)/s - 0.5L;
}

