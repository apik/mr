/* ****************************************************************
   rk6 (also rk6STU and rk6ST)
   ===========================

   Implements a 6-stage, 5th-order Cash-Karp-Fehlberg Runge-Kutta
   step.  The purpose of this (as opposed to classical 4-stage
   4th-order Runge-Kutta) is to enable efficient adaptive step size
   control.

   The criterion for deciding whether the step size is small enough is
   presently:

   |(5th order change) - (4th order embedded change)|
   < (precision_goal)*|5th order change|

   This is required for all dependent variables, unless the left-hand
   side is less than TSIL_TOL times the magnitude of the variable.
   Otherwise, the step is rejected, unless force_step == 1.

   Provided it finds the step size to be small enough, OR the argument
   force_step==1, then rk6 will return 1 and:

     if RKmode=0, then *RKindvar=s is incremented by an amount
     *RKdelta.

     if RKmode=1, then *RKindvar=ln(1-s/sthresh) is incremented by
     *RKdelta.

     If RKmode=2, then *RKindvar=ln(-s/qq) is incremented by *RKdelta.

   The dependent variables and their derivatives in *foo are also
   updated in preparation for the next step.

   If the step size needs to be decreased, then rk6 will return 0, and
   not increment the independent variable and not change the dependent
   variables or their derivatives.

   rk6 also replaces *RKdelta with an estimate of the optimal step
   size, regardless of whether that is going to be a retry of this
   step if it failed, or the next step if this one passed. But, it
   never tries to increase or decrease the step size by more than a
   factor of 2. And, the calling function Integrate() can and will
   reject the suggested step size if it gets to be too small or too
   large.

   ******************************************************************** */

#include "internal.h"

#define SafetyFactor 0.9L
#define TSIL_TINY    TSIL_TOL
#define TSIL_DIMLESS 1e-5
#define err_goal     0.01L*precision_goal /* This can be adjusted. */
#define max_allowed_error 1e4

/* Here are the Butcher coefficients for Cash-Karp-Fehlberg:*/

/* Confusingly, these are Numerical Recipes a_i */
#define ButchCKFc2 0.2L
#define ButchCKFc3 0.3L
#define ButchCKFc4 0.6L
#define ButchCKFc5 1.0L
#define ButchCKFc6 0.875L

/* Confusingly, these are Numerical Recipes b_ij */
#define ButchCKFa21 0.2L
#define ButchCKFa31 0.075L
#define ButchCKFa32 0.225L
#define ButchCKFa41 0.3L
#define ButchCKFa42 -0.9L
#define ButchCKFa43 1.2L
#define ButchCKFa51 -0.203703703703703703703703703704L
#define ButchCKFa52 2.5L
#define ButchCKFa53 -2.59259259259259259259259259259L
#define ButchCKFa54 1.29629629629629629629629629630L
#define ButchCKFa61 0.0294958043981481481481481481481L
#define ButchCKFa62 0.341796875L
#define ButchCKFa63 0.0415943287037037037037037037037L
#define ButchCKFa64 0.400345413773148148148148148148L
#define ButchCKFa65 0.061767578125L

/* Confusingly, these are Numerical Recipes c_i */
#define ButchCKFb1 0.0978835978835978835978835978836L
#define ButchCKFb2 0.0L
#define ButchCKFb3 0.402576489533011272141706924316L 
#define ButchCKFb4 0.210437710437710437710437710438L
#define ButchCKFb5 0.0L
#define ButchCKFb6 0.289102202145680406549971767363L

/* These are Numerical Recipes c_i - c_i^* */
#define ButchCKFe1 -0.00429377480158730158730158730159L
#define ButchCKFe2 0.0L 
#define ButchCKFe3 0.0186685860938578329882677708765L 
#define ButchCKFe4 -0.0341550268308080808080808080808L 
#define ButchCKFe5 -0.0193219866071428571428571428571L 
#define ButchCKFe6 0.0391022021456804065499717673631L 

/* **************************************************************** */

int rk6 (TSIL_DATA    *foo, 
         TSIL_COMPLEX *RKindvar,
         TSIL_COMPLEX *RKdelta,
         TSIL_REAL    sthresh,
         int          RKmode,
         TSIL_REAL    precision_goal, 
         int          force_step)
{
  static TSIL_COMPLEX startingB[2], startingS[2], startingT[6];
  static TSIL_COMPLEX startingU[4], startingM;
  TSIL_COMPLEX k1B[2], k2B[2], k3B[2], k4B[2], k5B[2], k6B[2];
  TSIL_COMPLEX k1S[2], k2S[2], k3S[2], k4S[2], k5S[2], k6S[2];
  TSIL_COMPLEX k1T[6], k2T[6], k3T[6], k4T[6], k5T[6], k6T[6];
  TSIL_COMPLEX k1U[4], k2U[4], k3U[4], k4U[4], k5U[4], k6U[4];
  TSIL_COMPLEX k1M,    k2M,    k3M,    k4M,    k5M,    k6M;
  TSIL_COMPLEX deltaB[2], deltaS[2], deltaT[6];
  TSIL_COMPLEX deltaU[4], deltaM;
  TSIL_REAL errorB[2], errorS[2], errorT[6];
  TSIL_REAL errorU[4], errorM;
  int i;
  TSIL_COMPLEX s, ds;
  TSIL_REAL qq = foo->qq;
  TSIL_COMPLEX next_step_size;
  TSIL_REAL temp1, temp2, maxerr; 
  TSIL_REAL step_rescale = 0.19L; /* This can only increase. */
  TSIL_REAL Sscale[2];
  int status = 0; 

  Sscale[0] = TSIL_CABS(foo->s) + (foo->v) + (foo->y) + (foo->z);
  Sscale[1] = TSIL_CABS(foo->s) + (foo->v) + (foo->x) + (foo->u);

  /* Relate s, ds to the independent variable. */
  if (0 == RKmode) {
    s = *RKindvar;
    ds = *RKdelta;
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar);
    ds = s * (*RKdelta);
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

  /* Fill up k1 arrays using existing values for derivatives: */
  for (i=0; i<2; i++) {
    k1B[i] = ds * (foo->B[i].deriv);
    k1S[i] = ds * (foo->S[i].deriv);
  }

  for (i=0; i<6; i++)
    k1T[i] = ds * (foo->T[i].deriv);

  for (i=0; i<4; i++)
    k1U[i] = ds * (foo->U[i].deriv);

  k1M = ds * (foo->M.deriv);

  /* BEGIN STAGE 2. */

  /* Update s, ds */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc2 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc2 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc2 * (*RKdelta) );
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + ButchCKFa21 * k1B[i];
    foo->S[i].value = startingS[i] + ButchCKFa21 * k1S[i];
  }

  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + ButchCKFa21 * k1T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + ButchCKFa21 * k1U[i];

  foo->M.value = startingM + ButchCKFa21 * k1M;

  /* Fill up k2 arrays: */
  for (i=0; i<2; i++) {
    k2B[i] = ds * dBds (foo->B[i], s);
    k2S[i] = ds * dSds (foo->S[i], s);
  }
  for (i=0; i<6; i++)
    k2T[i] = ds * dTds (foo->T[i], s);

  for (i=0; i<4; i++)
    k2U[i] = ds * dUds (foo->U[i], s);

  k2M = ds * dsMds (foo->M, s);

  /* BEGIN STAGE 3. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc3 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc3 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc3 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + ButchCKFa31 * k1B[i]
                                   + ButchCKFa32 * k2B[i];
    foo->S[i].value = startingS[i] + ButchCKFa31 * k1S[i]
                                   + ButchCKFa32 * k2S[i];
  }
  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + ButchCKFa31 * k1T[i]
                                   + ButchCKFa32 * k2T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + ButchCKFa31 * k1U[i]
                                   + ButchCKFa32 * k2U[i];

  foo->M.value = startingM + ButchCKFa31 * k1M
                           + ButchCKFa32 * k2M;

  /* Fill up k3 arrays: */
  for (i=0; i<2; i++) {
    k3B[i] = ds * dBds (foo->B[i] , s);
    k3S[i] = ds * dSds (foo->S[i] , s);
  }
  for (i=0; i<6; i++)
    k3T[i] = ds * dTds (foo->T[i] , s);

  for (i=0; i<4; i++)
    k3U[i] = ds * dUds (foo->U[i] , s);

  k3M = ds * dsMds (foo->M , s);

  /* BEGIN STAGE 4. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc4 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc4 * (*RKdelta)));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc4 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + ButchCKFa41 * k1B[i]
                                   + ButchCKFa42 * k2B[i]
                                   + ButchCKFa43 * k3B[i];
    foo->S[i].value = startingS[i] + ButchCKFa41 * k1S[i]
                                   + ButchCKFa42 * k2S[i]
                                   + ButchCKFa43 * k3S[i];
  }
  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + ButchCKFa41 * k1T[i]
                                   + ButchCKFa42 * k2T[i]
                                   + ButchCKFa43 * k3T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + ButchCKFa41 * k1U[i]
                                   + ButchCKFa42 * k2U[i]
                                   + ButchCKFa43 * k3U[i];

  foo->M.value = startingM + ButchCKFa41 * k1M
                           + ButchCKFa42 * k2M
                           + ButchCKFa43 * k3M;

  /* Fill up k4 arrays: */
  for (i=0; i<2; i++) {
    k4B[i] = ds * dBds (foo->B[i] , s);
    k4S[i] = ds * dSds (foo->S[i] , s);
  }
  for (i=0; i<6; i++)
    k4T[i] = ds * dTds (foo->T[i] , s);

  for (i=0; i<4; i++)
    k4U[i] = ds * dUds (foo->U[i] , s);

  k4M = ds * dsMds (foo->M , s);

  /* BEGIN STAGE 5. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc5 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc5 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc5 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + ButchCKFa51 * k1B[i]
                                   + ButchCKFa52 * k2B[i]
                                   + ButchCKFa53 * k3B[i]
                                   + ButchCKFa54 * k4B[i];
    foo->S[i].value = startingS[i] + ButchCKFa51 * k1S[i]
                                   + ButchCKFa52 * k2S[i]
                                   + ButchCKFa53 * k3S[i]
                                   + ButchCKFa54 * k4S[i];
  }
  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + ButchCKFa51 * k1T[i]
                                   + ButchCKFa52 * k2T[i]
                                   + ButchCKFa53 * k3T[i]
                                   + ButchCKFa54 * k4T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + ButchCKFa51 * k1U[i]
                                   + ButchCKFa52 * k2U[i]
                                   + ButchCKFa53 * k3U[i]
                                   + ButchCKFa54 * k4U[i];

  foo->M.value = startingM + ButchCKFa51 * k1M
                           + ButchCKFa52 * k2M
                           + ButchCKFa53 * k3M
                           + ButchCKFa54 * k4M;

  /* Fill up k5 arrays: */
  for (i=0; i<2; i++) {
    k5B[i] = ds * dBds (foo->B[i] , s);
    k5S[i] = ds * dSds (foo->S[i] , s);
  }
  for (i=0; i<6; i++)
    k5T[i] = ds * dTds (foo->T[i] , s);

  for (i=0; i<4; i++)
    k5U[i] = ds * dUds (foo->U[i] , s);

  k5M = ds * dsMds (foo->M , s);

  /* BEGIN STAGE 6. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc6 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc6 * (*RKdelta)));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc6 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  for (i=0; i<2; i++) {
    foo->B[i].value = startingB[i] + ButchCKFa61 * k1B[i]
                                   + ButchCKFa62 * k2B[i]
                                   + ButchCKFa63 * k3B[i]
                                   + ButchCKFa64 * k4B[i]
                                   + ButchCKFa65 * k5B[i];
    foo->S[i].value = startingS[i] + ButchCKFa61 * k1S[i]
                                   + ButchCKFa62 * k2S[i]
                                   + ButchCKFa63 * k3S[i]
                                   + ButchCKFa64 * k4S[i]
                                   + ButchCKFa65 * k5S[i];
  }

  for (i=0; i<6; i++)
    foo->T[i].value = startingT[i] + ButchCKFa61 * k1T[i]
                                   + ButchCKFa62 * k2T[i]
                                   + ButchCKFa63 * k3T[i]
                                   + ButchCKFa64 * k4T[i]
                                   + ButchCKFa65 * k5T[i];

  for (i=0; i<4; i++)
    foo->U[i].value = startingU[i] + ButchCKFa61 * k1U[i]
                                   + ButchCKFa62 * k2U[i]
                                   + ButchCKFa63 * k3U[i]
                                   + ButchCKFa64 * k4U[i]
                                   + ButchCKFa65 * k5U[i];

  foo->M.value = startingM + ButchCKFa61 * k1M
                           + ButchCKFa62 * k2M
                           + ButchCKFa63 * k3M
                           + ButchCKFa64 * k4M
                           + ButchCKFa65 * k5M;

  /* Fill up k6 arrays: */
  for (i=0; i<2; i++) {
    k6B[i] = ds * dBds (foo->B[i] , s);
    k6S[i] = ds * dSds (foo->S[i] , s);
  }
  for (i=0; i<6; i++)
    k6T[i] = ds * dTds (foo->T[i] , s);

  for (i=0; i<4; i++)
    k6U[i] = ds * dUds (foo->U[i] , s);

  k6M = ds * dsMds (foo->M , s);

  /* DONE WITH STAGES. */

  /* Compute 5th-order result changes in data values */
  for (i=0; i<2; i++) {
    deltaB[i] =   ButchCKFb1 * k1B[i] + ButchCKFb2 * k2B[i] 
                + ButchCKFb3 * k3B[i] + ButchCKFb4 * k4B[i] 
                + ButchCKFb5 * k5B[i] + ButchCKFb6 * k6B[i];
    deltaS[i] =   ButchCKFb1 * k1S[i] + ButchCKFb2 * k2S[i] 
                + ButchCKFb3 * k3S[i] + ButchCKFb4 * k4S[i] 
                + ButchCKFb5 * k5S[i] + ButchCKFb6 * k6S[i];
  }
  for (i=0; i<6; i++)
    deltaT[i] =   ButchCKFb1 * k1T[i] + ButchCKFb2 * k2T[i] 
                + ButchCKFb3 * k3T[i] + ButchCKFb4 * k4T[i] 
                + ButchCKFb5 * k5T[i] + ButchCKFb6 * k6T[i];
  for (i=0; i<4; i++)
    deltaU[i] =   ButchCKFb1 * k1U[i] + ButchCKFb2 * k2U[i] 
                + ButchCKFb3 * k3U[i] + ButchCKFb4 * k4U[i] 
                + ButchCKFb5 * k5U[i] + ButchCKFb6 * k6U[i];

  deltaM    =   ButchCKFb1 * k1M + ButchCKFb2 * k2M 
              + ButchCKFb3 * k3M + ButchCKFb4 * k4M
              + ButchCKFb5 * k5M + ButchCKFb6 * k6M;

  /* Compute estimated error of each dependent variable, and keep
     track of the maximum estimated step rescaling found. */

  maxerr = 0.0;

  for (i=0; i<2; i++) {
    errorB[i] = TSIL_CABS(ButchCKFe1 * k1B[i] + ButchCKFe2 * k2B[i] 
			+ ButchCKFe3 * k3B[i] + ButchCKFe4 * k4B[i] 
                        + ButchCKFe5 * k5B[i] + ButchCKFe6 * k6B[i]);

    temp1 = errorB[i]/(precision_goal*(TSIL_CABS(deltaB[i]) + TSIL_TINY));
    temp2 = errorB[i]/(err_goal*(TSIL_CABS(foo->B[i].value) + TSIL_DIMLESS));
    if (temp2 < temp1) temp1 = temp2;
    if (maxerr < temp1) maxerr = temp1;

    errorS[i] = TSIL_CABS(ButchCKFe1 * k1S[i] + ButchCKFe2 * k2S[i] 
                        + ButchCKFe3 * k3S[i] + ButchCKFe4 * k4S[i] 
                        + ButchCKFe5 * k5S[i] + ButchCKFe6 * k6S[i]);

    temp1 = errorS[i]/(precision_goal*(TSIL_CABS(deltaS[i]) + TSIL_TINY));
    temp2 = errorS[i]/(err_goal*(TSIL_CABS(foo->S[i].value) + TSIL_DIMLESS*Sscale[i]));
    if (temp2 < temp1) temp1 = temp2;
    if (maxerr < temp1) maxerr = temp1;
  }

  for (i=0; i<6; i++) {
    errorT[i] = TSIL_CABS(ButchCKFe1 * k1T[i] + ButchCKFe2 * k2T[i] 
                        + ButchCKFe3 * k3T[i] + ButchCKFe4 * k4T[i] 
                        + ButchCKFe5 * k5T[i] + ButchCKFe6 * k6T[i]);

    temp1 = errorT[i]/(precision_goal*(TSIL_CABS(deltaT[i]) + TSIL_TINY));
    temp2 = errorT[i]/(err_goal*(TSIL_CABS(foo->T[i].value) + TSIL_DIMLESS));
    if (temp2 < temp1) temp1 = temp2;
    if (maxerr < temp1) maxerr = temp1;
  }

  for (i=0; i<4; i++) {
    errorU[i] = TSIL_CABS(ButchCKFe1 * k1U[i] + ButchCKFe2 * k2U[i] 
                        + ButchCKFe3 * k3U[i] + ButchCKFe4 * k4U[i] 
                        + ButchCKFe5 * k5U[i] + ButchCKFe6 * k6U[i]);

    temp1 = errorU[i]/(precision_goal*(TSIL_CABS(deltaU[i]) + TSIL_TINY));
    temp2 = errorU[i]/(err_goal*(TSIL_CABS(foo->U[i].value) + TSIL_DIMLESS));
    if (temp2 < temp1) temp1 = temp2;
    if (maxerr < temp1) maxerr = temp1;
  }

  errorM = TSIL_CABS(ButchCKFe1 * k1M + ButchCKFe2 * k2M 
                   + ButchCKFe3 * k3M + ButchCKFe4 * k4M
                   + ButchCKFe5 * k5M + ButchCKFe6 * k6M);

    temp1 = errorM/(precision_goal*(TSIL_CABS(deltaM) + TSIL_TINY));
    temp2 = errorM/(err_goal*(TSIL_CABS(foo->M.value) + TSIL_DIMLESS));
    if (temp2 < temp1) temp1 = temp2;
    if (maxerr < temp1) maxerr = temp1;

  if (1 == TSIL_IsInfinite(deltaM)) maxerr = 10.L * max_allowed_error;

  for (i=0; i<4; i++) {
    if (1 == TSIL_IsInfinite(deltaU[i])) maxerr = 10.L * max_allowed_error;
  }

  for (i=0; i<6; i++) {
    if (1 == TSIL_IsInfinite(deltaT[i])) maxerr = 10.L * max_allowed_error;
  }

  for (i=0; i<2; i++) {
    if (1 == TSIL_IsInfinite(deltaS[i])) maxerr = 10.L * max_allowed_error;
  }

  for (i=0; i<2; i++) {
    if (1 == TSIL_IsInfinite(deltaB[i])) maxerr = 10.L * max_allowed_error;
  }

  /* If this step was forced, and the error was very large, and we're near 
     a threshold, then we're just not going to do any better so it is time 
     to bail out. Then the status returned is -1, which lets the calling 
     program Integrate know that it should stop immediately. */

  if ((1 == force_step) && (maxerr > max_allowed_error) && (1 == RKmode) &&
      (TSIL_CABS(*RKindvar) > 15))
  {
    status = -1;
  }

  /* Now, if the errors are acceptable, OR the step is being forced
     and the errors weren't huge, we increment the independent variable 
     and the data values and the derivatives, and get set to report 
     status = success. Otherwise, set the data back to their original 
     values.  
  */

  if ((maxerr < 1.0L) || ((1 == force_step) && (-1 != status))) 
  {      
    status = 1;

    /* Update data values */
    for (i=0; i<2; i++) {
      foo->B[i].value = startingB[i] + deltaB[i];
      foo->S[i].value = startingS[i] + deltaS[i]; 
    }
    for (i=0; i<6; i++)
      foo->T[i].value = startingT[i] + deltaT[i];
    for (i=0; i<4; i++)
      foo->U[i].value = startingU[i] + deltaU[i];
    foo->M.value = startingM + deltaM;

    /* Set independent variable up for next step. */
    *RKindvar += *RKdelta;

    /* Update independent variable. */
    if (0 == RKmode) {
      s = *RKindvar;
    }
    else if (1 == RKmode) {
      s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
    }
    else {
      s = -qq*TSIL_CEXP(*RKindvar);
    }

    /* Set derivatives for next step */
    for (i=0; i<2; i++) {
      foo->B[i].deriv = dBds (foo->B[i], s);
      foo->S[i].deriv = dSds (foo->S[i], s);
    }
    for (i=0; i<6; i++)
      foo->T[i].deriv = dTds (foo->T[i], s);

    for (i=0; i<4; i++)
      foo->U[i].deriv = dUds (foo->U[i], s);

    foo->M.deriv = dsMds (foo->M, s);
  }
  else {
    /* Set data back to original values, and don't touch the
       independent variable or the derivatives, since we failed to
       step. */
    for (i=0; i<2; i++) {
      foo->B[i].value =   startingB[i];
      foo->S[i].value =   startingS[i]; 
    }
    for (i=0; i<6; i++)
      foo->T[i].value =   startingT[i];
    for (i=0; i<4; i++)
      foo->U[i].value =   startingU[i];

    foo->M.value =   startingM;
  }

  /* Predict the appropriate next step size. */
  step_rescale = TSIL_SQRT(TSIL_SQRT(1.0L/maxerr));

  next_step_size = (*RKdelta) * (TSIL_COMPLEX) (SafetyFactor * step_rescale);

  /* Don't let the next step size be bigger than the present size by
     more than a factor of 1.5, or smaller than the present step size
     by more than a factor of 2. It isn't worth the risk. */

  if(TSIL_CABS(next_step_size) > 1.5L*TSIL_CABS(*RKdelta))
    next_step_size = 1.5L* (*RKdelta);

  if(TSIL_CABS(next_step_size) < 0.5L*TSIL_CABS(*RKdelta))
    next_step_size = 0.5L* (*RKdelta);

  /* Recommend the new step size to the calling function Integrate (). */
  *RKdelta = next_step_size;

  return status;
}

      
/* **************************************************************** */
/* This version takes a RK step for the STU case.                   */

int rk6_STU (TSIL_DATA   *foo, 
	    TSIL_COMPLEX *RKindvar,
	    TSIL_COMPLEX *RKdelta,
	    TSIL_REAL    sthresh,
	    int          RKmode,
	    TSIL_REAL    precision_goal, 
	    int          force_step)
{
  static TSIL_COMPLEX startingS, startingT[6];
  static TSIL_COMPLEX startingU;

  TSIL_COMPLEX k1S, k2S, k3S, k4S, k5S, k6S;
  /* Some wasted space here, but easier coding... */
  TSIL_COMPLEX k1T[6], k2T[6], k3T[6], k4T[6], k5T[6], k6T[6];
  TSIL_COMPLEX k1U, k2U, k3U, k4U, k5U, k6U;
  TSIL_COMPLEX deltaS, deltaT[6];
  TSIL_COMPLEX deltaU;
  TSIL_REAL errorS, errorT[6];
  TSIL_REAL errorU;
  int i;
  TSIL_COMPLEX s, ds;
  TSIL_REAL qq = foo->qq;
  TSIL_COMPLEX next_step_size;
  TSIL_REAL temp1, temp2, maxerr; 
  TSIL_REAL step_rescale = 0.19L; /* This can only increase. */
  TSIL_REAL Sscale;
  int status = 0; 

  /* Perhaps temporary: */
  int whichS = uxv;
  int whichU = xzuv;

/*   Sscale[0] = TSIL_CABS(foo->s) + (foo->v) + (foo->y) + (foo->z); */

  Sscale = TSIL_CABS(foo->s) + (foo->v) + (foo->x) + (foo->u);

  /* Relate s, ds to the independent variable. */
  if (0 == RKmode) {
    s = *RKindvar;
    ds = *RKdelta;
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar);
    ds = s * (*RKdelta);
  }

  /* Set the starting values */
  startingS = (foo->S[whichS].value);
  for (i=1; i<6; i+=2)
    startingT[i] = (foo->T[i].value);
  startingU = (foo->U[whichU].value);

  /* Fill up k1 arrays using existing values for derivatives: */
  k1S = ds * (foo->S[whichS].deriv);
  for (i=1; i<6; i+=2)
    k1T[i] = ds * (foo->T[i].deriv);
  k1U = ds * (foo->U[whichU].deriv);

  /* BEGIN STAGE 2. */

  /* Update s, ds */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc2 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc2 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc2 * (*RKdelta) );
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa21 * k1S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa21 * k1T[i];
  foo->U[whichU].value = startingU + ButchCKFa21 * k1U;

  /* Fill up k2 arrays: */
  k2S = ds * dSds (foo->S[whichS], s);
  for (i=1; i<6; i+=2)
    k2T[i] = ds * dTds (foo->T[i], s);
  k2U = ds * dUds (foo->U[whichU], s);

  /* BEGIN STAGE 3. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc3 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc3 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc3 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa31 * k1S
                                + ButchCKFa32 * k2S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa31 * k1T[i]
                                   + ButchCKFa32 * k2T[i];

  foo->U[whichU].value = startingU + ButchCKFa31 * k1U
                                 + ButchCKFa32 * k2U;

  /* Fill up k3 arrays: */
  k3S = ds * dSds (foo->S[whichS] , s);

  for (i=1; i<6; i+=2)
    k3T[i] = ds * dTds (foo->T[i] , s);

  k3U = ds * dUds (foo->U[whichU] , s);

  /* BEGIN STAGE 4. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc4 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc4 * (*RKdelta)));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc4 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa41 * k1S
                                + ButchCKFa42 * k2S
                                + ButchCKFa43 * k3S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa41 * k1T[i]
                                   + ButchCKFa42 * k2T[i]
                                   + ButchCKFa43 * k3T[i];

  foo->U[whichU].value = startingU + ButchCKFa41 * k1U
                                 + ButchCKFa42 * k2U
                                 + ButchCKFa43 * k3U;

  /* Fill up k4 arrays: */
  k4S = ds * dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k4T[i] = ds * dTds (foo->T[i] , s);
  k4U = ds * dUds (foo->U[whichU] , s);

  /* BEGIN STAGE 5. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc5 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc5 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc5 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa51 * k1S
                                + ButchCKFa52 * k2S
                                + ButchCKFa53 * k3S
                                + ButchCKFa54 * k4S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa51 * k1T[i]
                                   + ButchCKFa52 * k2T[i]
                                   + ButchCKFa53 * k3T[i]
                                   + ButchCKFa54 * k4T[i];

  foo->U[whichU].value = startingU + ButchCKFa51 * k1U
                                + ButchCKFa52 * k2U
                                + ButchCKFa53 * k3U
                                + ButchCKFa54 * k4U;

  /* Fill up k5 arrays: */
  k5S = ds * dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k5T[i] = ds * dTds (foo->T[i] , s);
  k5U = ds * dUds (foo->U[whichU] , s);

  /* BEGIN STAGE 6. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc6 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc6 * (*RKdelta)));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc6 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa61 * k1S
                                   + ButchCKFa62 * k2S
                                   + ButchCKFa63 * k3S
                                   + ButchCKFa64 * k4S
                                   + ButchCKFa65 * k5S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa61 * k1T[i]
                                   + ButchCKFa62 * k2T[i]
                                   + ButchCKFa63 * k3T[i]
                                   + ButchCKFa64 * k4T[i]
                                   + ButchCKFa65 * k5T[i];

  foo->U[whichU].value = startingU + ButchCKFa61 * k1U
                                   + ButchCKFa62 * k2U
                                   + ButchCKFa63 * k3U
                                   + ButchCKFa64 * k4U
                                   + ButchCKFa65 * k5U;

  /* Fill up k6 arrays: */
  k6S = ds * dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k6T[i] = ds * dTds (foo->T[i] , s);
  k6U = ds * dUds (foo->U[whichU] , s);

  /* DONE WITH STAGES. */

  /* Compute 5th-order result changes in data values */
  deltaS =   ButchCKFb1 * k1S + ButchCKFb2 * k2S 
           + ButchCKFb3 * k3S + ButchCKFb4 * k4S 
           + ButchCKFb5 * k5S + ButchCKFb6 * k6S;
  for (i=1; i<6; i+=2)
    deltaT[i] =   ButchCKFb1 * k1T[i] + ButchCKFb2 * k2T[i] 
                + ButchCKFb3 * k3T[i] + ButchCKFb4 * k4T[i] 
                + ButchCKFb5 * k5T[i] + ButchCKFb6 * k6T[i];
  deltaU =   ButchCKFb1 * k1U + ButchCKFb2 * k2U 
           + ButchCKFb3 * k3U + ButchCKFb4 * k4U 
           + ButchCKFb5 * k5U + ButchCKFb6 * k6U;

  /* Compute estimated error of each dependent variable, and keep
     track of the maximum estimated step rescaling found. */

  maxerr = 0.0;


  errorS = TSIL_CABS(ButchCKFe1 * k1S + ButchCKFe2 * k2S 
		     + ButchCKFe3 * k3S + ButchCKFe4 * k4S 
		     + ButchCKFe5 * k5S + ButchCKFe6 * k6S);

  temp1 = errorS/(precision_goal*(TSIL_CABS(deltaS) + TSIL_TINY));
  temp2 = errorS/(err_goal*(TSIL_CABS(foo->S[whichS].value) + TSIL_DIMLESS*Sscale));
  if (temp2 < temp1) temp1 = temp2;
  if (maxerr < temp1) maxerr = temp1;

  for (i=1; i<6; i+=2) {
    errorT[i] = TSIL_CABS(ButchCKFe1 * k1T[i] + ButchCKFe2 * k2T[i] 
                        + ButchCKFe3 * k3T[i] + ButchCKFe4 * k4T[i] 
                        + ButchCKFe5 * k5T[i] + ButchCKFe6 * k6T[i]);

    temp1 = errorT[i]/(precision_goal*(TSIL_CABS(deltaT[i]) + TSIL_TINY));
    temp2 = errorT[i]/(err_goal*(TSIL_CABS(foo->T[i].value) + TSIL_DIMLESS));
    if (temp2 < temp1) temp1 = temp2;
    if (maxerr < temp1) maxerr = temp1;
  }

  errorU = TSIL_CABS(ButchCKFe1 * k1U + ButchCKFe2 * k2U 
                        + ButchCKFe3 * k3U + ButchCKFe4 * k4U 
                        + ButchCKFe5 * k5U + ButchCKFe6 * k6U);

  temp1 = errorU/(precision_goal*(TSIL_CABS(deltaU) + TSIL_TINY));
  temp2 = errorU/(err_goal*(TSIL_CABS(foo->U[whichU].value) + TSIL_DIMLESS));
  if (temp2 < temp1) temp1 = temp2;
  if (maxerr < temp1) maxerr = temp1;

  if (1 == TSIL_IsInfinite(deltaU)) maxerr = 10.L * max_allowed_error;

  for (i=1; i<6; i+=2)
    if (1 == TSIL_IsInfinite(deltaT[i])) maxerr = 10.L * max_allowed_error;

  if (1 == TSIL_IsInfinite(deltaS)) maxerr = 10.L * max_allowed_error;

  /* If this step was forced, and the error was very large, and we're
     near a threshold, then we're just not going to do any better so
     it is time to bail out. Then the status returned is -1, which
     lets the calling program Integrate know that it should stop
     immediately. */

  if ((1 == force_step) && (maxerr > max_allowed_error) && (1 == RKmode) &&
      (TSIL_CABS(*RKindvar) > 15))
    {
      status = -1;
    }

  /* Now, if the errors are acceptable, OR the step is being forced
     and the errors weren't huge, we increment the independent
     variable and the data values and the derivatives, and get set to
     report status = success. Otherwise, set the data back to their
     original values.
  */

  if ((maxerr < 1.0L) || ((1 == force_step) && (-1 != status))) 
    {      
      status = 1;
      
      /* Update data values */
      foo->S[whichS].value = startingS + deltaS; 
      for (i=1; i<6; i+=2)
	foo->T[i].value = startingT[i] + deltaT[i];
      foo->U[whichU].value = startingU + deltaU;
      
      /* Set independent variable up for next step. */
      *RKindvar += *RKdelta;
      
      /* Update independent variable. */
      if (0 == RKmode) {
	s = *RKindvar;
      }
      else if (1 == RKmode) {
	s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
      }
      else {
	s = -qq*TSIL_CEXP(*RKindvar);
      }
      
      /* Set derivatives for next step */
      foo->S[whichS].deriv = dSds (foo->S[whichS], s);
      
      for (i=1; i<6; i+=2)
	foo->T[i].deriv = dTds (foo->T[i], s);
      
      foo->U[whichU].deriv = dUds (foo->U[whichU], s);
    }
  else {
    /* Set data back to original values, and don't touch the
       independent variable or the derivatives, since we failed to
       step. */
    foo->S[whichS].value = startingS; 
    
    for (i=1; i<6; i+=2)
      foo->T[i].value = startingT[i];
    foo->U[whichU].value = startingU;
  }

  /* Predict the appropriate next step size. */
  step_rescale = TSIL_SQRT(TSIL_SQRT(1.0L/maxerr));

  next_step_size = (*RKdelta) * (TSIL_COMPLEX) (SafetyFactor * step_rescale);

  /* Don't let the next step size be bigger than the present size by
     more than a factor of 1.5, or smaller than the present step size
     by more than a factor of 2. It isn't worth the risk. */

  if(TSIL_CABS(next_step_size) > 1.5L*TSIL_CABS(*RKdelta))
    next_step_size = 1.5L* (*RKdelta);

  if(TSIL_CABS(next_step_size) < 0.5L*TSIL_CABS(*RKdelta))
    next_step_size = 0.5L* (*RKdelta);

  /* Recommend the new step size to the calling function Integrate (). */
  *RKdelta = next_step_size;

  return status;
}


/* **************************************************************** */
/* This version takes a RK step for the ST case.                    */

int rk6_ST (TSIL_DATA *foo, 
	   TSIL_COMPLEX *RKindvar,
	   TSIL_COMPLEX *RKdelta,
	   TSIL_REAL sthresh,
	   int RKmode,
	   TSIL_REAL precision_goal, 
	   int force_step)
{
  static TSIL_COMPLEX startingS, startingT[6];
  TSIL_COMPLEX k1S, k2S, k3S, k4S, k5S, k6S;
  /* Some wasted space here, but easier coding... */
  TSIL_COMPLEX k1T[6], k2T[6], k3T[6], k4T[6], k5T[6], k6T[6];
  TSIL_COMPLEX deltaS, deltaT[6];
  TSIL_REAL errorS, errorT[6];
  int i;
  TSIL_COMPLEX s, ds;
  TSIL_REAL qq = foo->qq;
  TSIL_COMPLEX next_step_size;
  TSIL_REAL temp1, temp2, maxerr; 
  TSIL_REAL step_rescale = 0.19L; /* This can only increase. */
  TSIL_REAL Sscale;
  int status = 0; 

  /* Perhaps temporary: */
  int whichS = uxv;

/*   Sscale[0] = TSIL_CABS(foo->s) + (foo->v) + (foo->y) + (foo->z); */

  Sscale = TSIL_CABS(foo->s) + (foo->v) + (foo->x) + (foo->u);

  /* Relate s, ds to the independent variable. */
  if (0 == RKmode) {
    s = *RKindvar;
    ds = *RKdelta;
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar);
    ds = s * (*RKdelta);
  }

  /* Set the starting values */
  startingS = (foo->S[whichS].value);
  for (i=1; i<6; i+=2)
    startingT[i] = (foo->T[i].value);

  /* Fill up k1 arrays using existing values for derivatives: */
  k1S = ds * (foo->S[whichS].deriv);
  for (i=1; i<6; i+=2)
    k1T[i] = ds * (foo->T[i].deriv);

  /* BEGIN STAGE 2. */

  /* Update s, ds */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc2 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc2 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc2 * (*RKdelta) );
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa21 * k1S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa21 * k1T[i];

  /* Fill up k2 arrays: */
  k2S = ds * dSds (foo->S[whichS], s);
  for (i=1; i<6; i+=2)
    k2T[i] = ds * dTds (foo->T[i], s);

  /* BEGIN STAGE 3. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc3 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc3 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc3 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa31 * k1S
                                + ButchCKFa32 * k2S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa31 * k1T[i]
                                   + ButchCKFa32 * k2T[i];

  /* Fill up k3 arrays: */
  k3S = ds * dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k3T[i] = ds * dTds (foo->T[i] , s);

  /* BEGIN STAGE 4. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc4 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc4 * (*RKdelta)));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc4 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa41 * k1S
                                + ButchCKFa42 * k2S
                                + ButchCKFa43 * k3S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa41 * k1T[i]
                                   + ButchCKFa42 * k2T[i]
                                   + ButchCKFa43 * k3T[i];

  /* Fill up k4 arrays: */
  k4S = ds * dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k4T[i] = ds * dTds (foo->T[i] , s);

  /* BEGIN STAGE 5. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc5 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc5 * (*RKdelta) ) );
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc5 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa51 * k1S
                                + ButchCKFa52 * k2S
                                + ButchCKFa53 * k3S
                                + ButchCKFa54 * k4S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa51 * k1T[i]
                                   + ButchCKFa52 * k2T[i]
                                   + ButchCKFa53 * k3T[i]
                                   + ButchCKFa54 * k4T[i];

  /* Fill up k5 arrays: */
  k5S = ds * dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k5T[i] = ds * dTds (foo->T[i] , s);

  /* BEGIN STAGE 6. */

  /* Update independent variable. */
  if (0 == RKmode) {
    s = *RKindvar + ButchCKFc6 * (*RKdelta);
  }
  else if (1 == RKmode) {
    s = sthresh*(1.0L - TSIL_CEXP(*RKindvar + ButchCKFc6 * (*RKdelta)));
    ds = (s - sthresh) * (*RKdelta);
  }
  else {
    s = -qq*TSIL_CEXP(*RKindvar + ButchCKFc6 * (*RKdelta));
    ds = s * (*RKdelta);
  }

  /* Adjust data values */
  foo->S[whichS].value = startingS + ButchCKFa61 * k1S
                                   + ButchCKFa62 * k2S
                                   + ButchCKFa63 * k3S
                                   + ButchCKFa64 * k4S
                                   + ButchCKFa65 * k5S;
  for (i=1; i<6; i+=2)
    foo->T[i].value = startingT[i] + ButchCKFa61 * k1T[i]
                                   + ButchCKFa62 * k2T[i]
                                   + ButchCKFa63 * k3T[i]
                                   + ButchCKFa64 * k4T[i]
                                   + ButchCKFa65 * k5T[i];

  /* Fill up k6 arrays: */
  k6S = ds * dSds (foo->S[whichS] , s);
  for (i=1; i<6; i+=2)
    k6T[i] = ds * dTds (foo->T[i] , s);

  /* DONE WITH STAGES. */

  /* Compute 5th-order result changes in data values */
  deltaS =   ButchCKFb1 * k1S + ButchCKFb2 * k2S 
           + ButchCKFb3 * k3S + ButchCKFb4 * k4S 
           + ButchCKFb5 * k5S + ButchCKFb6 * k6S;

  for (i=1; i<6; i+=2)
    deltaT[i] =   ButchCKFb1 * k1T[i] + ButchCKFb2 * k2T[i] 
                + ButchCKFb3 * k3T[i] + ButchCKFb4 * k4T[i] 
                + ButchCKFb5 * k5T[i] + ButchCKFb6 * k6T[i];

  /* Compute estimated error of each dependent variable, and keep
     track of the maximum estimated step rescaling found. */

  maxerr = 0.0;

  errorS = TSIL_CABS(ButchCKFe1 * k1S + ButchCKFe2 * k2S 
         	   + ButchCKFe3 * k3S + ButchCKFe4 * k4S 
		   + ButchCKFe5 * k5S + ButchCKFe6 * k6S);

  temp1 = errorS/(precision_goal*(TSIL_CABS(deltaS) + TSIL_TINY));
  temp2 = errorS/(err_goal*(TSIL_CABS(foo->S[whichS].value) + TSIL_DIMLESS*Sscale));
  if (temp2 < temp1) temp1 = temp2;
  if (maxerr < temp1) maxerr = temp1;

  for (i=1; i<6; i+=2) {
    errorT[i] = TSIL_CABS(ButchCKFe1 * k1T[i] + ButchCKFe2 * k2T[i] 
                        + ButchCKFe3 * k3T[i] + ButchCKFe4 * k4T[i] 
                        + ButchCKFe5 * k5T[i] + ButchCKFe6 * k6T[i]);

    temp1 = errorT[i]/(precision_goal*(TSIL_CABS(deltaT[i]) + TSIL_TINY));
    temp2 = errorT[i]/(err_goal*(TSIL_CABS(foo->T[i].value) + TSIL_DIMLESS));
    if (temp2 < temp1) temp1 = temp2;
    if (maxerr < temp1) maxerr = temp1;
  }

  for (i=1; i<6; i+=2)
    if (1 == TSIL_IsInfinite(deltaT[i])) maxerr = 10.L * max_allowed_error;

  if (1 == TSIL_IsInfinite(deltaS)) maxerr = 10.L * max_allowed_error;

  /* If this step was forced, and the error was very large, and we're
     near a threshold, then we're just not going to do any better so
     it is time to bail out. Then the status returned is -1, which
     lets the calling program Integrate know that it should stop
     immediately. */

  if ((1 == force_step) && (maxerr > max_allowed_error) && (1 == RKmode) &&
      (TSIL_CABS(*RKindvar) > 15))
    {
      status = -1;
    }

  /* Now, if the errors are acceptable, OR the step is being forced
     and the errors weren't huge, we increment the independent
     variable and the data values and the derivatives, and get set to
     report status = success. Otherwise, set the data back to their
     original values.
  */

  if ((maxerr < 1.0L) || ((1 == force_step) && (-1 != status))) 
    {      
      status = 1;
      
      /* Update data values */
      foo->S[whichS].value = startingS + deltaS; 
      for (i=1; i<6; i+=2)
	foo->T[i].value = startingT[i] + deltaT[i];
      
      /* Set independent variable up for next step. */
      *RKindvar += *RKdelta;
      
      /* Update independent variable. */
      if (0 == RKmode) {
	s = *RKindvar;
      }
      else if (1 == RKmode) {
	s = sthresh*(1.0L - TSIL_CEXP(*RKindvar));
      }
      else {
	s = -qq*TSIL_CEXP(*RKindvar);
      }
      
      /* Set derivatives for next step */
      foo->S[whichS].deriv = dSds (foo->S[whichS], s);
      
      for (i=1; i<6; i+=2)
	foo->T[i].deriv = dTds (foo->T[i], s);
    }
  else {
    /* Set data back to original values, and don't touch the
       independent variable or the derivatives, since we failed to
       step. */
    foo->S[whichS].value = startingS; 
    
    for (i=1; i<6; i+=2)
      foo->T[i].value = startingT[i];
  }

  /* Predict the appropriate next step size. */
  step_rescale = TSIL_SQRT(TSIL_SQRT(1.0L/maxerr));

  next_step_size = (*RKdelta) * (TSIL_COMPLEX) (SafetyFactor * step_rescale);

  /* Don't let the next step size be bigger than the present size by
     more than a factor of 1.5, or smaller than the present step size
     by more than a factor of 2. It isn't worth the risk. */

  if(TSIL_CABS(next_step_size) > 1.5L*TSIL_CABS(*RKdelta))
    next_step_size = 1.5L* (*RKdelta);

  if(TSIL_CABS(next_step_size) < 0.5L*TSIL_CABS(*RKdelta))
    next_step_size = 0.5L* (*RKdelta);

  /* Recommend the new step size to the calling function Integrate (). */
  *RKdelta = next_step_size;

  return status;
}
