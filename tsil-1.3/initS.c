/* Routines for setting up and initializing S-type functions */

#include "internal.h"

/* ***************************************************************** */

void TSIL_ConstructS (TSIL_STYPE *s,
		 int n,
		 TSIL_REAL x,
		 TSIL_REAL y,
		 TSIL_REAL z, 
                 TSIL_REAL qq)
{
  s->which  = n;
  s->arg[0] = x;
  s->arg[1] = y;
  s->arg[2] = z;

  s->S_c = x + y + z - TSIL_A(x,qq) - TSIL_A(y,qq) - TSIL_A(z,qq);

  return;
}
