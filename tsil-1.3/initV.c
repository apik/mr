/* Routines for setup and initial values of V-type functions */

#include "internal.h"

/* **************************************************************** */

void TSIL_ConstructV (TSIL_VTYPE *V,
		 int n,
		 TSIL_REAL z,
		 TSIL_REAL x,
		 TSIL_REAL y,
		 TSIL_REAL v,
		 TSIL_REAL qq)
{
  V->which  = n;
  V->arg[0] = z;
  V->arg[1] = x;
  V->arg[2] = y;
  V->arg[3] = v;

  return;
}

