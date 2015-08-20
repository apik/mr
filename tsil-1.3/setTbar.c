/* Evaluation of Tbar functions. */

#include "internal.h"

/* ******************************************************************* */

void TSIL_SetTbar (TSIL_DATA *foo)
{
  TSIL_REAL x, y, z, u, v, qq;
  TSIL_COMPLEX s;

  if (foo->whichFns == STUM) {
    x = foo->x;
    y = foo->y;
    z = foo->z;
    u = foo->u;
    v = foo->v;
    s = foo->s;
    qq = foo->qq;
    
    if (v < TSIL_TOL) { 
      foo->T[vyz].value = TSIL_Infinity;
      foo->T[vxu].value = TSIL_Infinity;
      foo->Tbar[vyz].value = TSIL_Tbar0xy(y,z,s,qq);
      foo->Tbar[vxu].value = TSIL_Tbar0xy(x,u,s,qq);
    }
    else { 
      foo->Tbar[vyz].value = 
	(foo->T[vyz].value) + TSIL_LOG(v/qq)*TSIL_B(y,z,s,qq);
      foo->Tbar[vxu].value = 
	(foo->T[vxu].value) + TSIL_LOG(v/qq)*TSIL_B(x,u,s,qq);
    }
    
    if (x < TSIL_TOL) {
      foo->T[xuv].value = TSIL_Infinity;
      foo->Tbar[xuv].value = TSIL_Tbar0xy(u,v,s,qq);
    }
    else {
      foo->Tbar[xuv].value = 
	(foo->T[xuv].value) + TSIL_LOG(x/qq)*TSIL_B(u,v,s,qq);
    }
    
    if (y < TSIL_TOL) {
      foo->T[yzv].value = TSIL_Infinity;
      foo->Tbar[yzv].value = TSIL_Tbar0xy(z,v,s,qq);
    }
    else {
      foo->Tbar[yzv].value = 
	(foo->T[yzv].value) + TSIL_LOG(y/qq)*TSIL_B(z,v,s,qq);
    }
    
    if (z < TSIL_TOL) {
      foo->T[zyv].value = TSIL_Infinity;
      foo->Tbar[zyv].value = TSIL_Tbar0xy(y,v,s,qq);
    }
    else {
      foo->Tbar[zyv].value = 
	(foo->T[zyv].value) + TSIL_LOG(z/qq)*TSIL_B(y,v,s,qq);
    }
    
    if (u < TSIL_TOL) {
      foo->T[uxv].value = TSIL_Infinity;
      foo->Tbar[uxv].value = TSIL_Tbar0xy(x,v,s,qq);
    }
    else {
      foo->Tbar[uxv].value = 
	(foo->T[uxv].value) + TSIL_LOG(u/qq)*TSIL_B(x,v,s,qq);
    }
  }
  else {

    x = foo->x;
    u = foo->u;
    v = foo->v;
    s = foo->s;
    qq = foo->qq;
    
    if (v < TSIL_TOL) { 
      foo->T[vxu].value = TSIL_Infinity;
      foo->Tbar[vxu].value = TSIL_Tbar0xy(x,u,s,qq);
    }
    else { 
      foo->Tbar[vxu].value = 
	(foo->T[vxu].value) + TSIL_LOG(v/qq)*TSIL_B(x,u,s,qq);
    }
    
    if (x < TSIL_TOL) {
      foo->T[xuv].value = TSIL_Infinity;
      foo->Tbar[xuv].value = TSIL_Tbar0xy(u,v,s,qq);
    }
    else {
      foo->Tbar[xuv].value = 
	(foo->T[xuv].value) + TSIL_LOG(x/qq)*TSIL_B(u,v,s,qq);
    }
    
    if (u < TSIL_TOL) {
      foo->T[uxv].value = TSIL_Infinity;
      foo->Tbar[uxv].value = TSIL_Tbar0xy(x,v,s,qq);
    }
    else {
      foo->Tbar[uxv].value = 
	(foo->T[uxv].value) + TSIL_LOG(u/qq)*TSIL_B(x,v,s,qq);
    }
  }

  return;
}
