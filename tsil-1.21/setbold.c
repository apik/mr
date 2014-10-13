/* Evaluation of "bold" (UV divergent) forms of S,T,U functions. */

#include "internal.h"

/* **************************************************************** */

void SetBoldS (TSIL_STYPE *fun, TSIL_REAL s, TSIL_REAL qq)
{
  TSIL_REAL x, y, z;
  TSIL_REAL alphax, alphay, alphaz;

  /* For convenience */
  x = fun->arg[0];
  y = fun->arg[1];
  z = fun->arg[2];

  alphax = Alpha(x, qq);
  alphay = Alpha(y, qq);
  alphaz = Alpha(z, qq);

  /* 1/eps^0 term */
  fun->bold[0] = -0.5L*((x + y + z)*(1.0L + Zeta2) +
	       alphax*alphax + alphay*alphay + alphaz*alphaz) + fun->value;

  /* 1/eps^1 term */
  fun->bold[1] = 0.25L*s - 0.5L*(x + y + z) + A(x,qq) + A(y,qq) + A(z,qq);

  /* 1/eps^2 term */
  fun->bold[2] = -0.5*(x + y + z);

  return;
}

/* **************************************************************** */

void SetBoldT (TSIL_TTYPE *fun, TSIL_REAL qq)
{
  TSIL_REAL x, Axx;

  /* For convenience */
  x = fun->arg[0];

  if (x < TSIL_TOL) {
    fun->bold[0] = fun->bold[1] = TSIL_Infinity;
  }
  else {
    Axx = A(x, qq)/x;
    fun->bold[0] = 0.5L + Zeta2/2.0L + Axx + 0.5L*Axx*Axx + fun->value;
    fun->bold[1] = -0.5L - Axx;
  }

  fun->bold[2] = 0.5L;

  return;
}

/* **************************************************************** */

void SetBoldU (TSIL_UTYPE *fun, TSIL_REAL s, TSIL_REAL qq)
{
  TSIL_REAL x, y;

  /* DGR - added ss for type matching in calls below */
  TSIL_COMPLEX ss = s + 0.0L*I;

  /* For convenience */
  x = fun->arg[0];
  y = fun->arg[1];

/*   printf("\nSetting Bold U:\n"); */
/*   printf("U function = %Lf\n", fun->value); */
/*   printf("Beps       = ");TSIL_cprintf(Beps(x,y,ss,qq));printf("\n"); */

  fun->bold[0] = Beps(x, y, ss, qq) + fun->value;
  fun->bold[1] = 0.5L + B(x, y, ss, qq);
  fun->bold[2] = 0.5L;

  return;
}

/* **************************************************************** */

void SetBoldV (TSIL_VTYPE *fun, TSIL_REAL s, TSIL_REAL qq)
{
  TSIL_REAL x, y, Deltasxy;

  /* For convenience */
  x = fun->arg[0];
  y = fun->arg[1];
  Deltasxy = Delta(s,x,y);
  
  if (y/(x + y + TSIL_FABS(s)) < TSIL_TOL) {
    /* DGR commented out in v1.2 */
/*     TSIL_Warn("SetBoldV", "Vbold(x,y,z,u) is undefined for y = 0."); */
    fun->bold[0] = TSIL_Infinity;
    fun->bold[1] = TSIL_Infinity;
  }
  else if (TSIL_FABS((s - x - y - 2.0L*TSIL_SQRT(x*y))/(x + y)) < TSIL_TOL) {
    /* DGR commented out in v1.2 */
/*     TSIL_Warn("SetBoldV", "Vbold(x,y,z,u) is undefined for sqrt(s) = sqrt(x)+sqrt(y)."); */
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
    fun->bold[0] = ((s + x - y) * (Beps(x, y, s, qq) - 2.0L*B(x,y,s,qq))
                   + 2.0L * (Aeps(x, qq) - A(x,qq))
                   + (s - x - y) * (Aeps(y, qq) - A(y, qq))/y)/Deltasxy
                   + fun->value;
    fun->bold[1] = ((s + x - y) * (B(x,y,s,qq) - 1.0L) + 2.0L*A(x,qq)
                   + (s - x - y) * A(y,qq)/y)/Deltasxy;
  }

  fun->bold[2] = 0.0L;

  return;
}

/* **************************************************************** */

void SetBold (TSIL_DATA *foo)
{
  int i;

  if (foo->whichFns == STUM) {
    for (i=0; i<2; i++)
      SetBoldS (&foo->S[i], foo->s, foo->qq);

    for (i=0; i<6; i++)
      SetBoldT (&foo->T[i], foo->qq);

    for (i=0; i<4; i++)
      SetBoldU (&foo->U[i], foo->s, foo->qq);

    for (i=0; i<4; i++)
      SetBoldV (&foo->V[i], foo->s, foo->qq);
  }
  else if (foo->whichFns == STU) {

    SetBoldS (&foo->S[uxv], foo->s, foo->qq);

    for (i=1; i<6; i+=2)
      SetBoldT (&foo->T[i], foo->qq);

    SetBoldU (&foo->U[xzuv], foo->s, foo->qq);

    SetBoldV (&foo->V[xzuv], foo->s, foo->qq);
  }
  else if (foo->whichFns == ST) {

    SetBoldS (&foo->S[uxv], foo->s, foo->qq);

    for (i=1; i<6; i+=2)
      SetBoldT (&foo->T[i], foo->qq);
  }
  else
    TSIL_Error ("SetBold","This can't happen! whichFns not set",1);

  return;
}
