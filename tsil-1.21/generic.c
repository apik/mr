/* Evaluation for the generic case. */

#include "internal.h"
#include "tsil_params.h"
      
/* ******************************************************************* */
     
int MaxSteps (TSIL_DATA *foo, TSIL_COMPLEX z)
{
  return (foo->nStepsMaxCon 
          + floor((double) (TSIL_CABS(z) * foo->nStepsMaxVar)));
}

/* **************************************************************** */

void InitialValue (TSIL_DATA *foo, TSIL_COMPLEX sinit)
{
  TSIL_REAL x, y, z, u, v, qq;

  /* For convenience */
  x = foo->x;
  y = foo->y;
  z = foo->z;
  u = foo->u;
  v = foo->v;
  qq = foo->qq;

  /* DGR - Could probably reorganize this for greater clarity */
  if (foo->whichFns == STUM) {
    foo->B[xz].deriv = BprimeAtZero (x, z, qq);
    foo->B[xz].value = BAtZero (x, z, qq) + sinit * foo->B[xz].deriv;
    foo->B[yu].deriv = BprimeAtZero (y, u, qq);
    foo->B[yu].value = BAtZero (y, u, qq) + sinit * foo->B[yu].deriv;
  }

  if (foo->whichFns == STUM) {
    foo->S[vyz].deriv = SprimeAtZero (v, y, z, qq);
    foo->S[vyz].value = SAtZero (v, y, z, qq) + sinit * foo->S[vyz].deriv;
  }
  foo->S[uxv].deriv = SprimeAtZero (u, x, v, qq);
  foo->S[uxv].value = SAtZero (u, x, v, qq) + sinit * foo->S[uxv].deriv;

  /* We can set the T's to zero and disable their running if the first
     arg is zero, since they are TSIL_Infinity in this case, and don't
     enter into the running of other functions: */


  if (v > TSIL_TOL) {
    if (foo->whichFns == STUM) {
      foo->T[vyz].deriv = TprimeAtZero (v, y, z, qq);
      foo->T[vyz].value = TAtZero (v, y, z, qq) + sinit * foo->T[vyz].deriv;
    }
    foo->T[vxu].deriv = TprimeAtZero (v, x, u, qq);
    foo->T[vxu].value = TAtZero (v, x, u, qq) + sinit * foo->T[vxu].deriv;
  }
  else {
    foo->T[vyz].value = foo->T[vxu].value = 
      foo->T[vyz].deriv = foo->T[vxu].deriv = 0.0 + 0.0*I;
  }


  if (foo->whichFns == STUM) {
    if (y > TSIL_TOL) {
      foo->T[yzv].deriv = TprimeAtZero (y, z, v, qq);
      foo->T[yzv].value = TAtZero (y, z, v, qq) + sinit * foo->T[yzv].deriv;
    }
    else {
      foo->T[yzv].value = 
	foo->T[yzv].deriv = 0.0 + 0.0*I;
    }
  }

  if (foo->whichFns == STUM) {
    if (z > TSIL_TOL) {
      foo->T[zyv].deriv = TprimeAtZero (z, y, v, qq);
      foo->T[zyv].value = TAtZero (z, y, v, qq) + sinit * foo->T[zyv].deriv;
    }
    else {
      foo->T[zyv].value = 
	foo->T[zyv].deriv = 0.0 + 0.0*I;
    }
  }

  if (u > TSIL_TOL) {
    foo->T[uxv].deriv = TprimeAtZero (u, x, v, qq);
    foo->T[uxv].value = TAtZero (u, x, v, qq) + sinit * foo->T[uxv].deriv;
  }
  else {
    foo->T[uxv].value = 
      foo->T[uxv].deriv = 0.0 + 0.0*I;
  }

  if (x > TSIL_TOL) {
    foo->T[xuv].deriv = TprimeAtZero (x, u, v, qq);
    foo->T[xuv].value = TAtZero (x, u, v, qq) + sinit * foo->T[xuv].deriv;
  }
  else {
    foo->T[xuv].value = 
    foo->T[xuv].deriv = 0.0 + 0.0*I;
  }

  /* Similarly for the U's if the second arg is zero; here they are
     not infinite, but the running would otherwise give large roundoff
     errors and so they need correction at the end of the calculation
     anyway: */

  if (foo->whichFns == STUM) {
    if (x > TSIL_TOL) {
      foo->U[zxyv].deriv = UprimeAtZero (z, x, y, v, qq);
      foo->U[zxyv].value = UAtZero (z, x, y, v, qq) + sinit * foo->U[zxyv].deriv;
    }
    else {
      foo->U[zxyv].value = 
	foo->U[zxyv].deriv = 0.0 + 0.0*I;
    }

    if (y > TSIL_TOL) {
      foo->U[uyxv].deriv = UprimeAtZero (u, y, x, v, qq);
      foo->U[uyxv].value = UAtZero (u, y, x, v, qq) + sinit * foo->U[uyxv].deriv;
    }
    else {
      foo->U[uyxv].value = 
	foo->U[uyxv].deriv = 0.0 + 0.0*I;
    }
  }

  if (foo->whichFns != ST) {
    if (z > TSIL_TOL) {
      foo->U[xzuv].deriv = UprimeAtZero (x, z, u, v, qq);
      foo->U[xzuv].value = UAtZero (x, z, u, v, qq) + sinit * foo->U[xzuv].deriv;
    }
    else {
      foo->U[xzuv].value = 
	foo->U[xzuv].deriv = 0.0 + 0.0*I;
    }
  }

  if (foo->whichFns == STUM) {
    if (u > TSIL_TOL) {
      foo->U[yuzv].deriv = UprimeAtZero (y, u, z, v, qq);
      foo->U[yuzv].value = UAtZero (y, u, z, v, qq) + sinit * foo->U[yuzv].deriv;
    }
    else {
      foo->U[yuzv].value = 
	foo->U[yuzv].deriv = 0.0 + 0.0*I;
    }

    foo->M.deriv = sMprimeAtZero (x, y, z, u, v);
    foo->M.value = sMAtZero (x, y, z, u, v) + sinit * foo->M.deriv;
  }

  return;
}

/* **************************************************************** */
/* Designed to be used only for very small but non-zero sinit.      */

void InitialValueThreshAt0 (TSIL_DATA *foo, TSIL_COMPLEX sinit)
{
  TSIL_REAL x, y, z, u, v, qq;
  TSIL_COMPLEX rinit;

  /* For convenience */
  x = foo->x;
  y = foo->y;
  z = foo->z;
  u = foo->u;
  v = foo->v;
  qq = foo->qq;
  rinit = TSIL_CLOG(-sinit/qq);

  if (x + z < TSIL_TOL) {
    if (foo->whichFns == STUM) {
      foo->B[xz].value   = 2.0L - rinit;
      foo->B[xz].deriv   = -1.0L/sinit;    
    }
  }
  else {
    if (foo->whichFns == STUM) {
      foo->B[xz].value   = BAtZero (x, z, qq);
      foo->B[xz].deriv   = 0.0L;
      foo->U[zxyv].value = UAtZero (z, x, y, v, qq);
      foo->U[zxyv].deriv = 0.0L;
    }
    if (foo->whichFns != ST) {
      foo->U[xzuv].value = UAtZero (x, z, u, v, qq);
      foo->U[xzuv].deriv = 0.0L;
    }
  }

  if (foo->whichFns == STUM) {
    if (y + u < TSIL_TOL) {
      foo->B[yu].value   = 2.0L - rinit;
      foo->B[yu].deriv   = -1.0L/sinit;    
    }
    else { 
      foo->B[yu].value   = BAtZero (y, u, qq);
      foo->B[yu].deriv   = 0.0L;
      foo->U[yuzv].value = UAtZero (y, u, z, v, qq);
      foo->U[yuzv].deriv = 0.0L;
      foo->U[uyxv].value = UAtZero (u, y, x, v, qq);
      foo->U[uyxv].deriv = 0.0L;
    }
  }

  /* This is needed in all cases: */
  if (x + u + v < TSIL_TOL) {
    foo->S[uxv].value = sinit*(1.625L - 0.5L*rinit);
    foo->S[uxv].deriv = 1.125L - 0.5L*rinit;
  }
  else {
    foo->S[uxv].value = SAtZero (u, x, v, qq);
    foo->S[uxv].deriv = SprimeAtZero (u, x, v, qq);
    foo->T[uxv].value = TAtZero (u, x, v, qq);
    foo->T[xuv].value = TAtZero (x, u, v, qq);
    foo->T[vxu].value = TAtZero (v, x, u, qq);
  }

  if (foo->whichFns == STUM) {
    if (y + z + v < TSIL_TOL) {
      foo->S[vyz].value = sinit*(1.625L - 0.5L*rinit);
      foo->S[vyz].deriv = 1.125L - 0.5L*rinit;
    }
    else {
      foo->S[vyz].value = SAtZero (v, y, z, qq);
      foo->S[vyz].deriv = SprimeAtZero (v, y, z, qq);
      foo->T[vyz].value = TAtZero (v, y, z, qq);
      foo->T[yzv].value = TAtZero (y, z, v, qq);
      foo->T[zyv].value = TAtZero (z, y, v, qq);
    }
  }

  /* Next initialize M. */
  if (foo->whichFns == STUM)
    foo->M.value = 0.0L;

  /* Derivatives of T and M always initialized to zero when s=0 is a
     threshold, because dF/dr = s dF/ds -> 0 for very small s. */
  foo->T[uxv].deriv = 0.0L;
  foo->T[xuv].deriv = 0.0L;
  foo->T[vxu].deriv = 0.0L;
  if (foo->whichFns == STUM) {
    foo->T[vyz].deriv = 0.0L;
    foo->T[yzv].deriv = 0.0L;
    foo->T[zyv].deriv = 0.0L;
    foo->M.deriv = 0.0L;
  }

  /* We can set the T's to zero and disable their running if the first
     arg is zero, since they are TSIL_Infinity in this case.  (Note
     that the T derivatives are all taken care of above.)  Some of
     these are accomplished above so this is a bit redundant,
     but... */
  if (foo->whichFns == STUM) {
    if (v < TSIL_TOL)
      foo->T[vyz].value = foo->T[vxu].value = 0.0 + 0.0*I;

    if (y < TSIL_TOL)
      foo->T[yzv].value = 0.0 + 0.0*I;

    if (z < TSIL_TOL)
      foo->T[zyv].value = 0.0 + 0.0*I;
  }
  /* The rest are always needed: */
  if (u < TSIL_TOL)
    foo->T[uxv].value = 0.0 + 0.0*I;

  if (x < TSIL_TOL)
    foo->T[xuv].value = 0.0 + 0.0*I;

  /* Similarly for the U's if the second arg is zero; here they are
     not infinite, but the running gives garbage and they need
     correction at the end of the calculation anyway: */
  if (foo->whichFns == STUM) {
    if (x < TSIL_TOL) {
      foo->U[zxyv].value = 
	foo->U[zxyv].deriv = 0.0 + 0.0*I;
    }

    if (y < TSIL_TOL) {
      foo->U[uyxv].value = 
	foo->U[uyxv].deriv = 0.0 + 0.0*I;
    }

    if (u < TSIL_TOL) {
      foo->U[yuzv].value = 
	foo->U[yuzv].deriv = 0.0 + 0.0*I;
    }
  }
  if (foo->whichFns != ST) {
    if (z < TSIL_TOL) {
      foo->U[xzuv].value = 
	foo->U[xzuv].deriv = 0.0 + 0.0*I;
    }
  }

  return;
}

/* **************************************************************** */
/* Integration routine for STUM case.
    
   if intmode=0, integrate with independent variable s, assuming no
                 special points, ignore argument sthresh
    
   if intmode=1, integrate with independent variable s, assumes the
                 last point could be near a pseudo-threshold, not to
                 be used except for the last step of contour, ignore
                 argument sthresh
    
   if intmode=2, integrate with independent variable r = ln(1-s/sthresh),
                 assumes the last point could be near a threshold, not
                 to be used except for the last step of contour
  
   if intmode=3, integrate with independent variable r = ln(-s/qq),
                 assumes that s=0 could be a threshold,
                 ignores argument sthresh
    
   The minimum stepsize is     (t1 - t0)/max_steps, 
   The maximum stepsize is     (t1 - t0)/TSIL_nStepsMin, 
   The first guess stepsize is (t1 - t0)/TSIL_nStepsStart 
*/

int Integrate (TSIL_DATA    *foo,
               TSIL_COMPLEX t0,
               TSIL_COMPLEX t1,
               int          max_steps,
               int          intmode,
               TSIL_REAL    sthresh)
{
  TSIL_COMPLEX t, dt;
  TSIL_REAL pre_error = foo->precisionGoal;
  int force_step;
  int rkmode = 0;
  int goodsteps, badsteps; 
  int rk6status = 1; /* 1 for success or forced; 0 for need retry, 
                        -1 when the error is big and we need to bail. */

  goodsteps = badsteps = 0;

  if (intmode > 0)
    rkmode = intmode - 1;

  t  = t0;
  dt = (t1 - t0)/(foo->nStepsStart);
           
  /* Note in the following line, rk6status can be 1 or -1, but never 0 */    
  while ( (TSIL_CABS(dt) < 0.5*TSIL_CABS(t1-t)) && (1 == rk6status)) {
    for (;;) {
      if ( TSIL_CABS(dt) < TSIL_CABS(t1-t0)/max_steps ) {
        force_step = 1;
        dt = (t1-t0)/max_steps;
      }
      else force_step = 0;

      if (TSIL_CABS(dt) > TSIL_CABS(t1-t0)/(foo->nStepsMin))
        dt = (t1-t0)/(foo->nStepsMin);

      rk6status = foo->RKstepper6 (foo, &t, &dt, sthresh, rkmode, 
				   pre_error, force_step);

      if (0 != rk6status) { 
	goodsteps += (1-force_step);
	badsteps += force_step;
	break;
      }
    }
  }

  /* If the error got too big, it's over. Do not pass Go, do not
     collect $200. */
  if (-1 == rk6status) return 0;

  /* The remaining distance is less than twice the step size.  So, for
     the next-to-last step, go half the distance to the goal, and
     force it. Too many small steps here could be bad if there is a
     (pseudo)-threshold. */

  dt = 0.5L*(t1 - t);
  foo->RKstepper6 (foo, &t, &dt, sthresh, rkmode, pre_error, 1);

  /* Arrange final step to land exactly on t1, and force it. */
  dt = t1 - t;

  if ((0 == intmode) || (3 == intmode)) 
    foo->RKstepper6 (foo, &t, &dt, sthresh, rkmode, pre_error, 1);
  else 
    foo->RKstepper5 (foo, &t, dt, sthresh, rkmode);       

  /* Return a status code eventually */
  return 0;
}


/* **************************************************************** */
/* Handling of all generic (i.e. non-analytic) cases.               */

void CaseGeneric (TSIL_DATA *foo)
{
  TSIL_COMPLEX sInit, sFinal, rInit, rFinal, imDisp;
  TSIL_REAL    sthresh;
  TSIL_REAL    s = foo->s;
  TSIL_REAL    qq = foo->qq;
  TSIL_REAL    threshMin = foo->threshMin;
  TSIL_REAL    smallestspecialpoint;
  TSIL_REAL    temp;
  int          i;

  TSIL_Info("GENERIC CASE");

  /* Decide how to initialize; is s=0 a threshold, or close to one? */
  if (threshMin < TSIL_TOL) {
    TSIL_Info("There is a threshold at s=0.");
    sInit = I*SINIT;
    InitialValueThreshAt0 (foo, sInit);
  }
  else if (threshMin < THRESH_CUTOFF) {
    TSIL_Info("There is a threshold close to, but not at, s=0.");
    sInit = -SINIT;
    InitialValue (foo, sInit);
  }
  else {
    sInit = 0.L + 0.L*I;
    InitialValue (foo, 0.0L + 0.0L*I);
  }

  /* Find the point nearest s=0 that could give problems: */
  smallestspecialpoint = (foo->threshold)[0];

  for (i=1; i<(foo->nThresh); i++) {
    if ((foo->threshold)[i] < smallestspecialpoint) 
      smallestspecialpoint = (foo->threshold)[i];
  }

  for (i=0; i<(foo->nPthresh); i++) {
    if ((foo->pseudoThreshold)[i] < smallestspecialpoint) 
      smallestspecialpoint = (foo->pseudoThreshold)[i];
  }

  if (s < (smallestspecialpoint - THRESH_CUTOFF)) {
    /* Integrate along real s axis. */
    sFinal = (TSIL_COMPLEX) 0.5L*s;

    if (threshMin < THRESH_CUTOFF) {
      /* The smallest threshold is either 0 or close to 0, so change
         variables to r = lnbar(-s) for the first part of integration. */
      rInit  = TSIL_CLOG(-sInit/qq);
      temp = -0.5L*s/qq;
      if (temp > TSIL_TOL) rFinal = TSIL_CLOG(temp);
      else if (temp < -TSIL_TOL) rFinal = TSIL_CLOG(-temp) - I*PI;
      else rFinal = TSIL_CLOG(0.001L*TSIL_EPSILON) - I*0.5L*PI;
      Integrate (foo, rInit, rFinal, MaxSteps(foo,rFinal-rInit), 3, 0.0L);
    }
    else
      Integrate (foo, sInit, sFinal, MaxSteps(foo,sFinal-sInit), 0, 0.0L);

    sInit  = sFinal;
    sFinal = (TSIL_COMPLEX) s;
    Integrate (foo, sInit, sFinal, MaxSteps(foo, sFinal-sInit), 1, 0.0L);

    /* Set status value */
    foo->status = REAXIS;
  }
  else {
    /* Integrate in complex s plane.                            */
    /* No reason to go too far off the real axis if s is small. */
    if (s < IM_DISPL/10.0)
      imDisp = 10.0L * s * I;
    else
      imDisp = IM_DISPL * I;

    sFinal = imDisp;

    if (threshMin < THRESH_CUTOFF) {
      TSIL_Info("Using ln(-s/qq) as independent variable for first leg of contour.");
      rInit  = TSIL_CLOG(-sInit/qq);
      rFinal = TSIL_CLOG(-sFinal/qq);
      Integrate (foo, rInit, rFinal, MaxSteps(foo,rFinal-rInit), 3, 0.0L);
    }
    else Integrate (foo, sInit, sFinal, MaxSteps(foo,sFinal - sInit), 0, 0.0L);

    sInit  = sFinal;
    sFinal = s + imDisp;
    Integrate (foo, sInit, sFinal, MaxSteps(foo,sFinal - sInit), 0, 0.0L);

    sInit  = sFinal;
    sFinal = s;
    if (NearThreshold (foo, &sthresh, THRESH_CUTOFF) == YES) {
      if (TSIL_CABS(sthresh) < TSIL_TOL) {
        rInit  = TSIL_CLOG(-sInit/qq);
        rFinal = TSIL_CLOG(-sFinal/qq - I*TSIL_EPSILON);
        Integrate (foo, rInit, rFinal, MaxSteps(foo,rFinal-rInit), 3, 0.0L);
      }
      else {
        TSIL_Info("Using near-threshold stepper for final leg of contour.");
        rInit  = TSIL_CLOG(1.L - sInit/sthresh);
        temp = 1.L - s/sthresh;
        if (temp > TSIL_TOL) rFinal = TSIL_CLOG(temp);
        else if (temp < -TSIL_TOL) rFinal = TSIL_CLOG(-temp) - I*PI;
        else rFinal = TSIL_CLOG(0.001L*TSIL_EPSILON) - I*0.5L*PI;
        Integrate (foo, rInit, rFinal, MaxSteps(foo,rFinal - rInit), 2, sthresh);
      }
    }
    else 
      Integrate (foo, sInit, sFinal, MaxSteps(foo,sFinal - sInit), 1, 0.0L);

    /* Set status value */
    foo->status = CONTOUR;
  }

  /* Check whether we had a double pole case in any of the U's and fix
     it, if necessary: */
  if ((foo->x < TSIL_TOL) || (foo->y < TSIL_TOL) ||
      (foo->z < TSIL_TOL) || (foo->u < TSIL_TOL))
    CorrectUs (foo);

  /* Finally, convert s*M to M */
  foo->M.value /= s;

  return;
}
