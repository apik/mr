/* Analytic evaluations where possible. */

#include "internal.h"

/* ******************************************************************* */

int TSIL_CaseSpecial (TSIL_DATA *foo)
{
  TSIL_REAL x, y, z, u, v, s, qq;
  int success = 1;
  int tmpWarns;

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
  tmpWarns = printWarns;
  printWarns = NO;

  if (foo->whichFns == STUM) {
    foo->B[xz].value = TSIL_B(x, z, s, qq);
    foo->B[yu].value = TSIL_B(y, u, s, qq);

    success *= TSIL_Sanalytic (y, z, v, s, qq, &(foo->S[vyz].value));
    success *= TSIL_Sanalytic (u, x, v, s, qq, &(foo->S[uxv].value));
    
    success *= TSIL_Tanalytic (v, y, z, s, qq, &(foo->T[vyz].value));
    success *= TSIL_Tanalytic (v, x, u, s, qq, &(foo->T[vxu].value));
    success *= TSIL_Tanalytic (x, u, v, s, qq, &(foo->T[xuv].value));
    success *= TSIL_Tanalytic (y, z, v, s, qq, &(foo->T[yzv].value));
    success *= TSIL_Tanalytic (z, y, v, s, qq, &(foo->T[zyv].value));
    success *= TSIL_Tanalytic (u, x, v, s, qq, &(foo->T[uxv].value));
  
    success *= TSIL_Uanalytic (z, x, y, v, s, qq, &(foo->U[zxyv].value)); 
    success *= TSIL_Uanalytic (u, y, x, v, s, qq, &(foo->U[uyxv].value));
    success *= TSIL_Uanalytic (x, z, u, v, s, qq, &(foo->U[xzuv].value));
    success *= TSIL_Uanalytic (y, u, z, v, s, qq, &(foo->U[yuzv].value));
  }
  else if (foo->whichFns == STU) {
    foo->B[xz].value = TSIL_B(x, z, s, qq);

    success *= TSIL_Sanalytic (u, x, v, s, qq, &(foo->S[uxv].value));
    
    success *= TSIL_Tanalytic (v, x, u, s, qq, &(foo->T[vxu].value));
    success *= TSIL_Tanalytic (x, u, v, s, qq, &(foo->T[xuv].value));
    success *= TSIL_Tanalytic (u, x, v, s, qq, &(foo->T[uxv].value));
  
    success *= TSIL_Uanalytic (x, z, u, v, s, qq, &(foo->U[xzuv].value));
  }
  else {
    success *= TSIL_Sanalytic (u, x, v, s, qq, &(foo->S[uxv].value));

    success *= TSIL_Tanalytic (v, x, u, s, qq, &(foo->T[vxu].value));
    success *= TSIL_Tanalytic (x, u, v, s, qq, &(foo->T[xuv].value));
    success *= TSIL_Tanalytic (u, x, v, s, qq, &(foo->T[uxv].value));
  }

  /* Restore previous warnings setting: */
  printWarns = tmpWarns;

  return success;
}
