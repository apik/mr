/* Analytic evaluations where possible. */

#include "internal.h"

/* ******************************************************************* */

int CaseSpecial (TSIL_DATA *foo)
{
  TSIL_REAL x, y, z, u, v, s, qq;
  int success = 1;

  TSIL_Info("SPECIAL CASE");

  /* For convenience */
  x  = foo->x;
  y  = foo->y;
  z  = foo->z;
  u  = foo->u;
  v  = foo->v;
  s  = foo->s;
  qq = foo->qq;

  /* Temporarily disable WARNs */
  printWarns = NO;

  if (foo->whichFns == STUM) {
    foo->B[xz].value = B(x, z, s, qq);
    foo->B[yu].value = B(y, u, s, qq);

    success *= Sanalytic (y, z, v, s, qq, &(foo->S[vyz].value));
    success *= Sanalytic (u, x, v, s, qq, &(foo->S[uxv].value));
    
    success *= Tanalytic (v, y, z, s, qq, &(foo->T[vyz].value));
    success *= Tanalytic (v, x, u, s, qq, &(foo->T[vxu].value));
    success *= Tanalytic (x, u, v, s, qq, &(foo->T[xuv].value));
    success *= Tanalytic (y, z, v, s, qq, &(foo->T[yzv].value));
    success *= Tanalytic (z, y, v, s, qq, &(foo->T[zyv].value));
    success *= Tanalytic (u, x, v, s, qq, &(foo->T[uxv].value));
  
    success *= Uanalytic (z, x, y, v, s, qq, &(foo->U[zxyv].value)); 
    success *= Uanalytic (u, y, x, v, s, qq, &(foo->U[uyxv].value));
    success *= Uanalytic (x, z, u, v, s, qq, &(foo->U[xzuv].value));
    success *= Uanalytic (y, u, z, v, s, qq, &(foo->U[yuzv].value));
  }
  else if (foo->whichFns == STU) {
    foo->B[xz].value = B(x, z, s, qq);

    success *= Sanalytic (u, x, v, s, qq, &(foo->S[uxv].value));
    
    success *= Tanalytic (v, x, u, s, qq, &(foo->T[vxu].value));
    success *= Tanalytic (x, u, v, s, qq, &(foo->T[xuv].value));
    success *= Tanalytic (u, x, v, s, qq, &(foo->T[uxv].value));
  
    success *= Uanalytic (x, z, u, v, s, qq, &(foo->U[xzuv].value));
  }
  else {
    success *= Sanalytic (u, x, v, s, qq, &(foo->S[uxv].value));

    success *= Tanalytic (v, x, u, s, qq, &(foo->T[vxu].value));
    success *= Tanalytic (x, u, v, s, qq, &(foo->T[xuv].value));
    success *= Tanalytic (u, x, v, s, qq, &(foo->T[uxv].value));
  }

  /* Restore warnings */
  printWarns = YES;

  return success;
}
