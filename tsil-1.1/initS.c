/* Routines for setting up and initializing S-type functions */

#include "internal.h"

/* ***************************************************************** */

void ConstructS (TSIL_STYPE *s,
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

  s->S_c = x + y + z - A(x,qq) - A(y,qq) - A(z,qq);

  return;
}
