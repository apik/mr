/* Evaluation of "bold" (UV divergent) forms of S,T,U functions. */

#include "internal.h"

/* **************************************************************** */

void TSIL_SetBoldS (TSIL_STYPE *fun, TSIL_REAL s, TSIL_REAL qq)
{
  TSIL_REAL x, y, z;
  TSIL_REAL alphax, alphay, alphaz;

  /* For convenience */
  x = fun->arg[0];
  y = fun->arg[1];
  z = fun->arg[2];

  alphax = TSIL_Alpha(x, qq);
  alphay = TSIL_Alpha(y, qq);
  alphaz = TSIL_Alpha(z, qq);

  /* 1/eps^0 term */
  fun->bold[0] = -0.5L*((x + y + z)*(1.0L + Zeta2) +
	       alphax*alphax + alphay*alphay + alphaz*alphaz) + fun->value;

  /* 1/eps^1 term */
  fun->bold[1] = 0.25L*s - 0.5L*(x + y + z) + TSIL_A(x,qq) 
    + TSIL_A(y,qq) + TSIL_A(z,qq);

  /* 1/eps^2 term */
  fun->bold[2] = -0.5*(x + y + z);

  return;
}

/* **************************************************************** */

void TSIL_SetBoldT (TSIL_TTYPE *fun, TSIL_REAL qq)
{
  TSIL_REAL x, Axx;

  /* For convenience */
  x = fun->arg[0];

  if (x < TSIL_TOL) {
    fun->bold[0] = fun->bold[1] = TSIL_Infinity;
  }
  else {
    Axx = TSIL_A(x, qq)/x;
    fun->bold[0] = 0.5L + Zeta2/2.0L + Axx + 0.5L*Axx*Axx + fun->value;
    fun->bold[1] = -0.5L - Axx;
  }

  fun->bold[2] = 0.5L;

  return;
}

/* **************************************************************** */

void TSIL_SetBoldU (TSIL_UTYPE *fun, TSIL_REAL s, TSIL_REAL qq)
{
  TSIL_REAL x, y;

  /* DGR - added ss for type matching in calls below */
  TSIL_COMPLEX ss = s + 0.0L*I;

  /* For convenience */
  x = fun->arg[0];
  y = fun->arg[1];

  fun->bold[0] = TSIL_Beps(x, y, ss, qq) + fun->value;
  fun->bold[1] = 0.5L + TSIL_B(x, y, ss, qq);
  fun->bold[2] = 0.5L;

  return;
}

/* **************************************************************** */

void TSIL_SetBoldV (TSIL_VTYPE *fun, TSIL_REAL s, TSIL_REAL qq)
{
  TSIL_REAL x, y, Deltasxy;

  /* For convenience */
  x = fun->arg[0];
  y = fun->arg[1];
  Deltasxy = TSIL_Delta(s,x,y);
  
  if (y/(x + y + TSIL_FABS(s)) < TSIL_TOL) {
    /* DGR commented out in v1.2 */
/*     TSIL_Warn("TSIL_SetBoldV", "Vbold(x,y,z,u) is undefined for y = 0."); */
    fun->bold[0] = TSIL_Infinity;
    fun->bold[1] = TSIL_Infinity;
  }
  else if (TSIL_FABS((s - x - y - 2.0L*TSIL_SQRT(x*y))/(x + y)) < TSIL_TOL) {
    /* DGR commented out in v1.2 */
/*     TSIL_Warn("TSIL_SetBoldV", "Vbold(x,y,z,u) is undefined for sqrt(s) = sqrt(x)+sqrt(y)."); */
    fun->bold[0] = TSIL_Infinity;
    fun->bold[1] = TSIL_Infinity;
  }
  else if (TSIL_FABS((s - x - y + 2.0L*TSIL_SQRT(x*y))/(x + y)) < TSIL_TOL) {
    if (TSIL_FABS(x - y)/(x + y) > TSIL_TOL) {
      fun->bold[0] = fun->value + (2.0L + TSIL_CLOG(x/qq) + 
                     (-2.0L - TSIL_CLOG(y/qq))*TSIL_CSQRT(x/y) +
                     0.25L*(TSIL_CLOG(x/qq)*TSIL_CLOG(x/qq) - 
                            TSIL_CLOG(y/qq)*TSIL_CLOG(y/qq)))/s;
      fun->bold[1] = (TSIL_SQRT(x/y) - 1.0L + 0.5L*TSIL_LOG(y/x))/
                     (x + y - 2.0L*TSIL_SQRT(x*y));
    }
    else {
      fun->bold[0] = fun->value - 0.5L*TSIL_CLOG(x/qq)/x;
      fun->bold[1] = 0.5L/x;
    }
  }
  else {
    fun->bold[0] = ((s + x - y) * (TSIL_Beps(x, y, s, qq)
				   - 2.0L*TSIL_B(x,y,s,qq))
               + 2.0L * (TSIL_Aeps(x, qq) - TSIL_A(x,qq))
               + (s - x - y) * (TSIL_Aeps(y, qq) - TSIL_A(y, qq))/y)/Deltasxy
                   + fun->value;
    fun->bold[1] = ((s + x - y) * (TSIL_B(x,y,s,qq) - 1.0L) + 2.0L*TSIL_A(x,qq)
                   + (s - x - y) * TSIL_A(y,qq)/y)/Deltasxy;
  }

  fun->bold[2] = 0.0L;

  return;
}

/* **************************************************************** */

void TSIL_SetBold (TSIL_DATA *foo)
{
  int i;

  if (foo->whichFns == STUM) {
    for (i=0; i<2; i++)
      TSIL_SetBoldS (&foo->S[i], foo->s, foo->qq);

    for (i=0; i<6; i++)
      TSIL_SetBoldT (&foo->T[i], foo->qq);

    for (i=0; i<4; i++)
      TSIL_SetBoldU (&foo->U[i], foo->s, foo->qq);

    for (i=0; i<4; i++)
      TSIL_SetBoldV (&foo->V[i], foo->s, foo->qq);
  }
  else if (foo->whichFns == STU) {

    TSIL_SetBoldS (&foo->S[uxv], foo->s, foo->qq);

    for (i=1; i<6; i+=2)
      TSIL_SetBoldT (&foo->T[i], foo->qq);

    TSIL_SetBoldU (&foo->U[xzuv], foo->s, foo->qq);

    TSIL_SetBoldV (&foo->V[xzuv], foo->s, foo->qq);
  }
  else if (foo->whichFns == ST) {

    TSIL_SetBoldS (&foo->S[uxv], foo->s, foo->qq);

    for (i=1; i<6; i+=2)
      TSIL_SetBoldT (&foo->T[i], foo->qq);
  }
  else
    TSIL_Error ("TSIL_SetBold","This can't happen! whichFns not set",1);

  return;
}
