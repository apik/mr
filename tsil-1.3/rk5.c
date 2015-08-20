/* 
  Implements a general 5-stage, 4th-order Runge-Kutta step.  The
  purpose of this (as opposed to classical 4-substep 4th order
  Runge-Kutta) is to never have to evaluate anything at the endpoint.
  Butcher coefficients a[i,j], b[i], c[i] have to be selected
  according to theory; for our purposes it is important that c[5] < 1.

  This code can be made more efficient by removing terms in the code
  involving coefficient a[i,j] that are 0. However, efficiency is not
  a big concern, because the function TSIL_rk5 should only be called once
  in each run. So all of the 0 coefficients are left in as a nod to
  generality.

  If RKmode = 0, then RKindvar = s is incremented by an amount
  RKdelta.  (Here sthresh is ignored.)

  If RKmode = 1, then RKindvar = ln(1-s/sthresh) is incremented by
  RKdelta.  (Here sthresh is the nearest threshold.)

  The component members of *foo and *RKindvar are updated at the end.
  The derivatives are NOT updated, because they might diverge causing
  a major speed hit, and because they are not needed (since this is
  ONLY to be used for the last step!)

  NOTE: do NOT call rk4 AFTER calling TSIL_rk5, because rk4 expects the
  derivatives to be up-to-date in the struct, and TSIL_rk5 doesn't update
  them!  It is permitted to call TSIL_rk5 repeatedly as a test, because TSIL_rk5
  doesn't assume that the derivatives in the struct are updated. But
  that would be needlessly inefficient.  */

#include "internal.h"

/* Here is a nice set of Butcher coefficients. Others are possible.*/
#define Butcherc2 0.25L
#define Butcherc3 0.375L
#define Butcherc4 0.5L
#define Butcherc5 0.625L

#define Butchera21 0.25L
#define Butchera31 0.0L
#define Butchera32 0.375L
#define Butchera41 0.0L
#define Butchera42 0.5L
#define Butchera43 0.0L
#define Butchera51 0.0L
#define Butchera52 35.0L/72.0L
#define Butchera53 0.0L
#define Butchera54 (5.0L/36.0L)

#define Butcherb1 (-1.0L/15.0L)
#define Butcherb2 (2.0L/3.0L)
#define Butcherb3 (4.0L/3.0L)
#define Butcherb4 (-10.0L/3.0L)
#define Butcherb5 (12.0L/5.0L)

/* **************************************************************** */

void TSIL_rk5 (TSIL_DATA *foo, 
	       TSIL_COMPLEX *RKindvar,
	       TSIL_COMPLEX RKdelta,
	       TSIL_REAL sthresh,
	       int RKmode)
{
  static TSIL_COMPLEX k1B[2], k2B[2], k3B[2], k4B[2], k5B[2];
  static TSIL_COMPLEX k1S[2], k2S[2], k3S[2], k4S[2], k5S[2];
  static TSIL_COMPLEX k1T[6], k2T[6], k3T[6], k4T[6], k5T[6];
  static TSIL_COMPLEX k1U[4], k2U[4], k3U[4], k4U[4], k5U[4];
  static TSIL_COMPLEX k1M,    k2M,    k3M,    k4M,    k5M;
  static TSIL_COMPLEX startingB[2], startingS[2], startingT[6];
  static TSIL_COMPLEX startingU[4], startingM;

  int i;
  TSIL_COMPLEX s, ds;

  /* Relate s, ds to the independent variable. */
  if (0 == RKmode) {
    s = *RKindvar;
    ds = RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
    ds = (s - sthresh)* RKdelta;
  }

  /* Set the starting values */
  for (i=0; i<2; i++) {
    startingB[i] = (foo->B[i].value);
    startingS[i] = (foo->S[i].value);
  }
  for (i=0; i<6; i++) 
    startingT[i] = (foo->T[i].value);

  for (i=0; i<4; i++)
    startingU[i] = (foo->U[i].value);

  startingM = (foo->M.value);

  /* Fill up k1 arrays: */
  for (i=0; i<2; i++) {
    k1B[i] = ds * TSIL_dBds_rk (foo->B[i], s);
    k1S[i] = ds * TSIL_dSds (foo->S[i], s);
  }
  for (i=0; i<6; i++)
    k1T[i] = ds * TSIL_dTds (foo->T[i], s);

  for (i=0; i<4; i++)
    k1U[i] = ds * TSIL_dUds (foo->U[i], s);

  k1M = ds * TSIL_dsMds (foo->M, s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc2 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc2 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + Butchera21 * k1B[i];
    foo->S[i].value = startingS[i] + Butchera21 * k1S[i];
  }
  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + Butchera21 * k1T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + Butchera21 * k1U[i];

  foo->M.value = startingM + Butchera21 * k1M;

  /* Fill up k2 arrays: */
  for (i=0; i<2; i++) {
    k2B[i] = ds * TSIL_dBds_rk (foo->B[i], s);
    k2S[i] = ds * TSIL_dSds (foo->S[i], s);
  }
  for (i=0; i<6; i++)
    k2T[i] = ds * TSIL_dTds (foo->T[i], s);

  for (i=0; i<4; i++)
    k2U[i] = ds * TSIL_dUds (foo->U[i], s);

  k2M = ds * TSIL_dsMds (foo->M, s);

  /* Update independent variabl. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc3 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc3 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + Butchera31 * k1B[i]
                                   + Butchera32 * k2B[i];
    foo->S[i].value = startingS[i] + Butchera31 * k1S[i]
                                   + Butchera32 * k2S[i];
  }
  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + Butchera31 * k1T[i]
                                   + Butchera32 * k2T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + Butchera31 * k1U[i]
                                   + Butchera32 * k2U[i];

  foo->M.value = startingM + Butchera31 * k1M
                           + Butchera32 * k2M;

  /* Fill up k3 arrays: */
  for (i=0; i<2; i++) {
    k3B[i] = ds * TSIL_dBds_rk (foo->B[i] , s);
    k3S[i] = ds * TSIL_dSds (foo->S[i] , s);
  }
  for (i=0; i<6; i++)
    k3T[i] = ds * TSIL_dTds (foo->T[i] , s);

  for (i=0; i<4; i++)
    k3U[i] = ds * TSIL_dUds (foo->U[i] , s);

  k3M = ds * TSIL_dsMds (foo->M , s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc4 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc4 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + Butchera41 * k1B[i]
                                   + Butchera42 * k2B[i]
                                   + Butchera43 * k3B[i];
    foo->S[i].value = startingS[i] + Butchera41 * k1S[i]
                                   + Butchera42 * k2S[i]
                                   + Butchera43 * k3S[i];
  }
  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + Butchera41 * k1T[i]
                                   + Butchera42 * k2T[i]
                                   + Butchera43 * k3T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + Butchera41 * k1U[i]
                                   + Butchera42 * k2U[i]
                                   + Butchera43 * k3U[i];

  foo->M.value = startingM + Butchera41 * k1M
                           + Butchera42 * k2M
                           + Butchera43 * k3M;

  /* Fill up k4 arrays: */
  for (i=0; i<2; i++) {
    k4B[i] = ds * TSIL_dBds_rk (foo->B[i] , s);
    k4S[i] = ds * TSIL_dSds (foo->S[i] , s);
  }
  for (i=0; i<6; i++)
    k4T[i] = ds * TSIL_dTds (foo->T[i] , s);

  for (i=0; i<4; i++)
    k4U[i] = ds * TSIL_dUds (foo->U[i] , s);

  k4M = ds * TSIL_dsMds (foo->M , s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc5 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc5 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + Butchera51 * k1B[i]
                                   + Butchera52 * k2B[i]
                                   + Butchera53 * k3B[i]
                                   + Butchera54 * k4B[i];
    foo->S[i].value = startingS[i] + Butchera51 * k1S[i]
                                   + Butchera52 * k2S[i]
                                   + Butchera53 * k3S[i]
                                   + Butchera54 * k4S[i];
  }
  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + Butchera51 * k1T[i]
                                   + Butchera52 * k2T[i]
                                   + Butchera53 * k3T[i]
                                   + Butchera54 * k4T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + Butchera51 * k1U[i]
                                   + Butchera52 * k2U[i]
                                   + Butchera53 * k3U[i]
                                   + Butchera54 * k4U[i];

  foo->M.value = startingM + Butchera51 * k1M
                           + Butchera52 * k2M
                           + Butchera53 * k3M
                           + Butchera54 * k4M;

  /* Fill up k5 arrays: */
  for (i=0; i<2; i++) {
    k5B[i] = ds * TSIL_dBds_rk (foo->B[i] , s);
    k5S[i] = ds * TSIL_dSds (foo->S[i] , s);
  }
  for (i=0; i<6; i++)
    k5T[i] = ds * TSIL_dTds (foo->T[i] , s);

  for (i=0; i<4; i++)
    k5U[i] = ds * TSIL_dUds (foo->U[i] , s);

  k5M = ds * TSIL_dsMds (foo->M , s);

  /* Increment data values */
  for (i=0; i<2; i++) {
    foo->B[i].value =   startingB[i]
                      + Butcherb1 * k1B[i] + Butcherb2 * k2B[i] 
                      + Butcherb3 * k3B[i] + Butcherb4 * k4B[i] 
                      + Butcherb5 * k5B[i];
    foo->S[i].value =   startingS[i] 
                      + Butcherb1 * k1S[i] + Butcherb2 * k2S[i] 
                      + Butcherb3 * k3S[i] + Butcherb4 * k4S[i] 
                      + Butcherb5 * k5S[i];
  }
  for (i=0; i<6; i++)
    foo->T[i].value =   startingT[i]
                      + Butcherb1 * k1T[i] + Butcherb2 * k2T[i] 
                      + Butcherb3 * k3T[i] + Butcherb4 * k4T[i] 
                      + Butcherb5 * k5T[i];
  for (i=0; i<4; i++)
    foo->U[i].value =   startingU[i]
                      + Butcherb1 * k1U[i] + Butcherb2 * k2U[i] 
                      + Butcherb3 * k3U[i] + Butcherb4 * k4U[i] 
                      + Butcherb5 * k5U[i];

  foo->M.value =   startingM
                 + Butcherb1 * k1M + Butcherb2 * k2M 
                 + Butcherb3 * k3M + Butcherb4 * k4M
                 + Butcherb5 * k5M;

  /* Update independent variable for next step, if there is one. */
  *RKindvar += RKdelta;

  return;
}


/* **************************************************************** */
/* This versions take a RK step for the STU case.                   */

void TSIL_rk5_STU (TSIL_DATA    *foo, 
		   TSIL_COMPLEX *RKindvar,
		   TSIL_COMPLEX RKdelta,
		   TSIL_REAL    sthresh,
		   int          RKmode)
{
  static TSIL_COMPLEX k1S, k2S, k3S, k4S, k5S;
  static TSIL_COMPLEX k1T[6], k2T[6], k3T[6], k4T[6], k5T[6];
  static TSIL_COMPLEX k1U, k2U, k3U, k4U, k5U;

  static TSIL_COMPLEX startingS, startingT[6];
  static TSIL_COMPLEX startingU;

  int i;
  TSIL_COMPLEX s, ds;

  int whichS = uxv;
  int whichU = xzuv;

  /* Relate s, ds to the independent variable. */
  if (0 == RKmode) {
    s = *RKindvar;
    ds = RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
    ds = (s - sthresh)* RKdelta;
  }

  /* Set the starting values */
  startingS = (foo->S[whichS].value);
  for (i=1; i<6; i+=2) 
    startingT[i] = (foo->T[i].value);

  startingU = (foo->U[whichU].value);

  /* Fill up k1 arrays: */
  k1S = ds * TSIL_dSds (foo->S[whichS], s);
  for (i=1; i<6; i+=2)
    k1T[i] = ds * TSIL_dTds (foo->T[i], s);
  k1U = ds * TSIL_dUds (foo->U[whichU], s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc2 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc2 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + Butchera21 * k1S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + Butchera21 * k1T[i];
  foo->U[whichU].value = startingU + Butchera21 * k1U;

  /* Fill up k2 arrays: */
  k2S = ds * TSIL_dSds (foo->S[whichS], s);
  for (i=1; i<6; i+=2)
    k2T[i] = ds * TSIL_dTds (foo->T[i], s);
  k2U = ds * TSIL_dUds (foo->U[whichU], s);

  /* Update independent variabl. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc3 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc3 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + Butchera31 * k1S
                                   + Butchera32 * k2S;

  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + Butchera31 * k1T[i]
                                   + Butchera32 * k2T[i];

  foo->U[whichU].value = startingU + Butchera31 * k1U
                                   + Butchera32 * k2U;

  /* Fill up k3 arrays: */
  k3S = ds * TSIL_dSds (foo->S[whichS] , s);

  for (i=1; i<6; i+=2)
    k3T[i] = ds * TSIL_dTds (foo->T[i] , s);

  k3U = ds * TSIL_dUds (foo->U[whichU] , s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc4 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc4 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + Butchera41 * k1S
                                   + Butchera42 * k2S
                                   + Butchera43 * k3S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + Butchera41 * k1T[i]
                                   + Butchera42 * k2T[i]
                                   + Butchera43 * k3T[i];

  foo->U[whichU].value = startingU + Butchera41 * k1U
                                   + Butchera42 * k2U
                                   + Butchera43 * k3U;

  /* Fill up k4 arrays: */
  k4S = ds * TSIL_dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k4T[i] = ds * TSIL_dTds (foo->T[i] , s);
  k4U = ds * TSIL_dUds (foo->U[whichU] , s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc5 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc5 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + Butchera51 * k1S
                                   + Butchera52 * k2S
                                   + Butchera53 * k3S
                                   + Butchera54 * k4S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + Butchera51 * k1T[i]
                                   + Butchera52 * k2T[i]
                                   + Butchera53 * k3T[i]
                                   + Butchera54 * k4T[i];

  foo->U[whichU].value = startingU + Butchera51 * k1U
                                   + Butchera52 * k2U
                                   + Butchera53 * k3U
                                   + Butchera54 * k4U;

  /* Fill up k5 arrays: */
  k5S = ds * TSIL_dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k5T[i] = ds * TSIL_dTds (foo->T[i] , s);
  k5U = ds * TSIL_dUds (foo->U[whichU] , s);

  /* Increment data values */
  foo->S[whichS].value =   startingS
                      + Butcherb1 * k1S + Butcherb2 * k2S 
                      + Butcherb3 * k3S + Butcherb4 * k4S 
                      + Butcherb5 * k5S;

  for (i=1; i<6; i+=2)
    foo->T[i].value =   startingT[i]
                      + Butcherb1 * k1T[i] + Butcherb2 * k2T[i] 
                      + Butcherb3 * k3T[i] + Butcherb4 * k4T[i] 
                      + Butcherb5 * k5T[i];
  
  foo->U[whichU].value =   startingU
                      + Butcherb1 * k1U + Butcherb2 * k2U 
                      + Butcherb3 * k3U + Butcherb4 * k4U 
                      + Butcherb5 * k5U;

  /* Update independent variable for next step, if there is one. */
  *RKindvar += RKdelta;

  return;
}


/* **************************************************************** */
/* This versions take a RK step for the ST case.                    */

void TSIL_rk5_ST (TSIL_DATA *foo, 
		  TSIL_COMPLEX *RKindvar,
		  TSIL_COMPLEX RKdelta,
		  TSIL_REAL sthresh,
		  int RKmode)
{
  static TSIL_COMPLEX k1S, k2S, k3S, k4S, k5S;
  static TSIL_COMPLEX k1T[6], k2T[6], k3T[6], k4T[6], k5T[6];
  static TSIL_COMPLEX startingS, startingT[6];
  int i;
  TSIL_COMPLEX s, ds;
  int whichS = uxv;

  /* Relate s, ds to the independent variable. */
  if (0 == RKmode) {
    s = *RKindvar;
    ds = RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
    ds = (s - sthresh)* RKdelta;
  }

  /* Set the starting values */
  startingS = (foo->S[whichS].value);
  for (i=1; i<6; i+=2) 
    startingT[i] = (foo->T[i].value);

  /* Fill up k1 arrays: */
  k1S = ds * TSIL_dSds (foo->S[whichS], s);
  for (i=1; i<6; i+=2)
    k1T[i] = ds * TSIL_dTds (foo->T[i], s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc2 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc2 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + Butchera21 * k1S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + Butchera21 * k1T[i];

  /* Fill up k2 arrays: */
  k2S = ds * TSIL_dSds (foo->S[whichS], s);
  for (i=1; i<6; i+=2)
    k2T[i] = ds * TSIL_dTds (foo->T[i], s);

  /* Update independent variabl. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc3 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc3 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + Butchera31 * k1S
                                   + Butchera32 * k2S;

  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + Butchera31 * k1T[i]
                                   + Butchera32 * k2T[i];

  /* Fill up k3 arrays: */
  k3S = ds * TSIL_dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k3T[i] = ds * TSIL_dTds (foo->T[i] , s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc4 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc4 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + Butchera41 * k1S
                                   + Butchera42 * k2S
                                   + Butchera43 * k3S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + Butchera41 * k1T[i]
                                   + Butchera42 * k2T[i]
                                   + Butchera43 * k3T[i];

  /* Fill up k4 arrays: */
  k4S = ds * TSIL_dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k4T[i] = ds * TSIL_dTds (foo->T[i] , s);

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + Butcherc5 * RKdelta;
  }
  else {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + Butcherc5 * RKdelta));
    ds = (s - sthresh) * RKdelta;
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + Butchera51 * k1S
                                   + Butchera52 * k2S
                                   + Butchera53 * k3S
                                   + Butchera54 * k4S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + Butchera51 * k1T[i]
                                   + Butchera52 * k2T[i]
                                   + Butchera53 * k3T[i]
                                   + Butchera54 * k4T[i];

  /* Fill up k5 arrays: */
  k5S = ds * TSIL_dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k5T[i] = ds * TSIL_dTds (foo->T[i] , s);

  /* Increment data values */
  foo->S[whichS].value =   startingS
                      + Butcherb1 * k1S + Butcherb2 * k2S 
                      + Butcherb3 * k3S + Butcherb4 * k4S 
                      + Butcherb5 * k5S;

  for (i=1; i<6; i+=2)
    foo->T[i].value =   startingT[i]
                      + Butcherb1 * k1T[i] + Butcherb2 * k2T[i] 
                      + Butcherb3 * k3T[i] + Butcherb4 * k4T[i] 
                      + Butcherb5 * k5T[i];
  
  /* Update independent variable for next step, if there is one. */
  *RKindvar += RKdelta;

  return;
}
