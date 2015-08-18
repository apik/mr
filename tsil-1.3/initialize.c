/*
  Initialization (pointer alignment) in general data struct.

  This is everything that is independent of the values of the
  arguments; everything that *does* depend on x,y,z,u,v,qq is set in
  TSIL_SetParameters below.
*/

#include "internal.h"
#include "tsil_params.h"
int printWarns = YES;
/* FILE *warnfile, *errfile; */

/* **************************************************************** */

void TSIL_Construct (TSIL_DATA *foo)
{
  int i, j;

  /* These define which functions go where in the evolution equations: */
  static int SrefT[2][3] = {{Tvyz, Tyzv, Tzyv},
			    {Tuxv, Txuv, Tvxu}};

  static int TrefS[6]    = {Svyz, Suxv, Svyz, Suxv, Svyz, Suxv};
  static int TrefT[6][3] = {{Tvyz, Tyzv, Tzyv},
			    {Tuxv, Txuv, Tvxu},
			    {Tyzv, Tzyv, Tvyz},
			    {Txuv, Tuxv, Tvxu},
			    {Tzyv, Tyzv, Tvyz},
			    {Tvxu, Txuv, Tuxv}};

  static int UrefS[4]    = {Svyz, Suxv, Suxv, Svyz};
  static int UrefT[4][3] = {{Tzyv, Tvyz, Tyzv},
			    {Tuxv, Tvxu, Txuv},
			    {Txuv, Tvxu, Tuxv},
			    {Tyzv, Tvyz, Tzyv}};

  static int MrefB[2] = {Bxz, Byu};
  static int MrefS[2] = {Suxv, Svyz};
  static int MrefT[6] = {Txuv, Tyzv, Tzyv, Tuxv, Tvxu, Tvyz};
  static int MrefU[4] = {Uxzuv, Uyuzv, Uzxyv, Uuyxv};

  /* Align pointers: */

  /* S-type objects */
  for (i=0; i<2; i++)   /* Loop over S objects */
    for (j=0; j<3; j++) /* Loop over pointers in S */
      foo->S[i].tval[j] = &(foo->T[SrefT[i][j]].value);

  /* T-type objects */
  for (i=0; i<6; i++) { /* Loop over T objects */
    foo->T[i].sval = &(foo->S[TrefS[i]].value);
    for (j=0; j<3; j++) /* Loop over pointers in T */
      foo->T[i].tval[j] = &(foo->T[TrefT[i][j]].value);
  }

  /* U-type objects */
  for (i=0; i<4; i++) { /* Loop over U objects */
    foo->U[i].sval = &(foo->S[UrefS[i]].value);
    for (j=0; j<3; j++) /* Loop over pointers in U */
      foo->U[i].tval[j] = &(foo->T[UrefT[i][j]].value);
  }

  /* M-type object */
  for (i=0; i<2; i++) {
    foo->M.bval[i] = &(foo->B[MrefB[i]].value);
    foo->M.sval[i] = &(foo->S[MrefS[i]].value);
  }
  for (i=0; i<6; i++)
    foo->M.tval[i] = &(foo->T[MrefT[i]].value);

  for (i=0; i<4; i++)
    foo->M.uval[i] = &(foo->U[MrefU[i]].value);

  foo->isAligned = YES;

  return;
}

/* **************************************************************** */
/* Sets parameters (x,y,z,u,v,qq) in a data struct and updates all  */
/* sub-objects and evolution coefficients.                          */

int TSIL_SetParameters (TSIL_DATA *foo,
			TSIL_REAL x,
			TSIL_REAL y,
			TSIL_REAL z,
			TSIL_REAL u,
			TSIL_REAL v,
			TSIL_REAL qq)
{
  TSIL_REAL tmp[6];
  int i, mtype;

  /* Basic sanity check on arguments: */
  tmp[0] = x;
  tmp[1] = y;
  tmp[2] = z;
  tmp[3] = u;
  tmp[4] = v;
  tmp[5] = qq;

  for (i=0; i<5; i++)
    if (tmp[i] < 0.0)
      TSIL_Error("TSIL_SetParameters",
		 "Squared mass argument cannot be negative", 2);

  if (tmp[5] < TSIL_TOL*TSIL_TOL)
    TSIL_Error("TSIL_SetParameters",
	       "Renormalization scale squared must be positive", 2);

  /* Set up data object if necessary */
  /* if (foo->isAligned != YES) */
  TSIL_Construct (foo);
  /* printWarns = YES; */

  /* STUM evaluation */
  foo->whichFns = STUM;
  foo->RKstepper6 = &TSIL_rk6;
  foo->RKstepper5 = &TSIL_rk5;

  /* Set values in the data object */
  foo->x  = x;
  foo->y  = y;
  foo->z  = z;
  foo->u  = u;
  foo->v  = v;
  foo->qq = qq;

  /* Set default scale factor */
  foo->scaleFac = 1.0L;

  /* Calculate thresholds... */
  foo->threshold[0] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(z), 2);
  foo->threshold[1] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(u), 2);
  foo->threshold[2] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(u) + TSIL_SQRT(v), 2);
  foo->threshold[3] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(z) + TSIL_SQRT(v), 2);
  foo->threshMin    = TSIL_MinAbs (foo->threshold, 4);
  foo->nThresh = 4;

  /* ...and pseudo-thresholds: */
  foo->pseudoThreshold[0] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(z), 2);
  foo->pseudoThreshold[1] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(u), 2);
  foo->pseudoThreshold[2] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(u) - TSIL_SQRT(v), 2);
  foo->pseudoThreshold[3] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(u) + TSIL_SQRT(v), 2);
  foo->pseudoThreshold[4] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(u) - TSIL_SQRT(v), 2);
  foo->pseudoThreshold[5] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(z) - TSIL_SQRT(v), 2);
  foo->pseudoThreshold[6] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(z) + TSIL_SQRT(v), 2);
  foo->pseudoThreshold[7] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(z) - TSIL_SQRT(v), 2);
  foo->psThreshMin        = TSIL_MinAbs (foo->pseudoThreshold, 8);
  foo->nPthresh = 8;

  /* Set values in all sub-objects: */
  TSIL_ConstructB (&(foo->B[xz]), x, z, qq);
  TSIL_ConstructB (&(foo->B[yu]), y, u, qq);

  TSIL_ConstructS (&(foo->S[vyz]), vyz, v, y, z, qq);
  TSIL_ConstructS (&(foo->S[uxv]), uxv, u, x, v, qq);

  TSIL_ConstructT (&(foo->T[vyz]), vyz, v, y, z, qq);
  TSIL_ConstructT (&(foo->T[uxv]), uxv, u, x, v, qq);
  TSIL_ConstructT (&(foo->T[yzv]), yzv, y, z, v, qq);
  TSIL_ConstructT (&(foo->T[xuv]), xuv, x, u, v, qq);
  TSIL_ConstructT (&(foo->T[zyv]), zyv, z, y, v, qq);
  TSIL_ConstructT (&(foo->T[vxu]), vxu, v, x, u, qq);

  TSIL_ConstructU (&(foo->U[zxyv]), zxyv, z, x, y, v, qq);
  TSIL_ConstructU (&(foo->U[uyxv]), uyxv, u, y, x, v, qq);
  TSIL_ConstructU (&(foo->U[xzuv]), xzuv, x, z, u, v, qq);
  TSIL_ConstructU (&(foo->U[yuzv]), yuzv, y, u, z, v, qq);

  TSIL_ConstructV (&(foo->V[zxyv]), zxyv, z, x, y, v, qq);
  TSIL_ConstructV (&(foo->V[uyxv]), uyxv, u, y, x, v, qq);
  TSIL_ConstructV (&(foo->V[xzuv]), xzuv, x, z, u, v, qq);
  TSIL_ConstructV (&(foo->V[yuzv]), yuzv, y, u, z, v, qq);

  mtype = TSIL_ConstructM (&(foo->M), x, y, z, u, v, qq);

  /* Set Runge-Kutta step-size parameter parameters.  */
  /* These can be reset by calling TSIL_ResetStepSizeParams */

  foo->precisionGoal = TSIL_PRECISION_GOAL;
  foo->nStepsStart   = TSIL_NSTEPS_START;
  foo->nStepsMin     = TSIL_NSTEPS_MIN;
  foo->nStepsMaxCon  = TSIL_NSTEPS_MAX_CON;
  foo->nStepsMaxVar  = TSIL_NSTEPS_MAX_VAR;

  foo->isInitialized = TRUE;
  foo->status        = UNEVALUATED;

  return 0;
}

/* **************************************************************** */
/* Sets parameters (x,y,z,u,v,qq) for the case STU and updates all  */
/* sub-objects and evolution coefficients.                          */

int TSIL_SetParametersSTU (TSIL_DATA *foo,
			   TSIL_REAL x,
			   TSIL_REAL z,
			   TSIL_REAL u,
			   TSIL_REAL v,
			   TSIL_REAL qq)
{
  TSIL_REAL tmp[5];
  int i;

  /* Basic sanity check on arguments: */
  tmp[0] = x;
  tmp[1] = z;
  tmp[2] = u;
  tmp[3] = v;
  tmp[4] = qq;

  for (i=0; i<4; i++)
    if (tmp[i] < 0.0)
      TSIL_Error("TSIL_SetParametersSTU",
		 "Squared mass argument cannot be negative", 2);

  if (tmp[4] < TSIL_TOL*TSIL_TOL)
    TSIL_Error("TSIL_SetParametersSTU",
	       "Renormalization scale squared must be positive", 2);

  /* Set up data object if necessary */
  /* if (foo->isAligned != YES) */
  TSIL_Construct (foo);

  foo->whichFns = STU;
  foo->RKstepper6 = &TSIL_rk6_STU;
  foo->RKstepper5 = &TSIL_rk5_STU;

  /* Set values in the data object */
  /* DGR -- leave y undefined, or set to nan or something? */
  foo->x  = x;
  foo->y  = TSIL_Infinity;
  foo->z  = z;
  foo->u  = u;
  foo->v  = v;
  foo->qq = qq;

  /* Set default scale factor */
  foo->scaleFac = 1.0L;

  /* Calculate thresholds... */
  foo->threshold[0] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(z), 2);
/*   foo->threshold[1] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(u), 2); */
  foo->threshold[1] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(u) + TSIL_SQRT(v), 2);
/*   foo->threshold[3] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(z) + TSIL_SQRT(v), 2); */
  foo->threshMin    = TSIL_MinAbs (foo->threshold, 2);
  foo->nThresh = 2;

  /* ...and pseudo-thresholds: */
  foo->pseudoThreshold[0] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(z), 2);
/*   foo->pseudoThreshold[1] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(u), 2); */
  foo->pseudoThreshold[1] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(u) - TSIL_SQRT(v), 2);
  foo->pseudoThreshold[2] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(u) + TSIL_SQRT(v), 2);
  foo->pseudoThreshold[3] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(u) - TSIL_SQRT(v), 2);
/*   foo->pseudoThreshold[5] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(z) - TSIL_SQRT(v), 2); */
/*   foo->pseudoThreshold[6] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(z) + TSIL_SQRT(v), 2); */
/*   foo->pseudoThreshold[7] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(z) - TSIL_SQRT(v), 2); */
  foo->psThreshMin        = TSIL_MinAbs (foo->pseudoThreshold, 4);
  foo->nPthresh = 4;

  /* Set values in all sub-objects: */
/*   TSIL_ConstructB (&(foo->B[xz]), x, z, qq); */
/*   TSIL_ConstructB (&(foo->B[yu]), y, u, qq); */

/*   TSIL_ConstructS (&(foo->S[vyz]), vyz, v, y, z, qq); */
  TSIL_ConstructS (&(foo->S[uxv]), uxv, u, x, v, qq);

/*   TSIL_ConstructT (&(foo->T[vyz]), vyz, v, y, z, qq); */
  TSIL_ConstructT (&(foo->T[uxv]), uxv, u, x, v, qq);
/*   TSIL_ConstructT (&(foo->T[yzv]), yzv, y, z, v, qq); */
  TSIL_ConstructT (&(foo->T[xuv]), xuv, x, u, v, qq);
/*   TSIL_ConstructT (&(foo->T[zyv]), zyv, z, y, v, qq); */
  TSIL_ConstructT (&(foo->T[vxu]), vxu, v, x, u, qq);

/*   TSIL_ConstructU (&(foo->U[zxyv]), zxyv, z, x, y, v, qq); */
/*   TSIL_ConstructU (&(foo->U[uyxv]), uyxv, u, y, x, v, qq); */
  TSIL_ConstructU (&(foo->U[xzuv]), xzuv, x, z, u, v, qq);
/*   TSIL_ConstructU (&(foo->U[yuzv]), yuzv, y, u, z, v, qq); */

/*   TSIL_ConstructV (&(foo->V[zxyv]), zxyv, z, x, y, v, qq); */
/*   TSIL_ConstructV (&(foo->V[uyxv]), uyxv, u, y, x, v, qq); */
  TSIL_ConstructV (&(foo->V[xzuv]), xzuv, x, z, u, v, qq);
/*   TSIL_ConstructV (&(foo->V[yuzv]), yuzv, y, u, z, v, qq); */

/*   mtype = TSIL_ConstructM (&(foo->M), x, y, z, u, v, qq); */

  /* Set Runge-Kutta step-size parameter parameters.  */
  /* These can be reset by calling TSIL_ResetStepSizeParams */

  foo->precisionGoal = TSIL_PRECISION_GOAL;
  foo->nStepsStart   = TSIL_NSTEPS_START;
  foo->nStepsMin     = TSIL_NSTEPS_MIN;
  foo->nStepsMaxCon  = TSIL_NSTEPS_MAX_CON;
  foo->nStepsMaxVar  = TSIL_NSTEPS_MAX_VAR;

  foo->isInitialized = TRUE;
  foo->status        = UNEVALUATED;

  return 0;
}

/* **************************************************************** */
/* Sets parameters (x,y,z,u,v,qq) for the case ST and updates all   */
/* sub-objects and evolution coefficients.                          */

int TSIL_SetParametersST (TSIL_DATA *foo,
			  TSIL_REAL x,
			  TSIL_REAL u,
			  TSIL_REAL v,
			  TSIL_REAL qq)
{
  TSIL_REAL tmp[4];
  int i;

  /* Basic sanity check on arguments: */
  tmp[0] = x;
  tmp[1] = u;
  tmp[2] = v;
  tmp[3] = qq;

  for (i=0; i<3; i++)
    if (tmp[i] < 0.0)
      TSIL_Error("TSIL_SetParametersST",
		 "Squared mass argument cannot be negative", 2);

  if (tmp[3] < TSIL_TOL*TSIL_TOL)
    TSIL_Error("TSIL_SetParametersST",
	       "Renormalization scale squared must be positive", 2);

  /* Set up data object if necessary */
  /* if (foo->isAligned != YES) */
  TSIL_Construct (foo);

  foo->whichFns = ST;
  foo->RKstepper6 = &TSIL_rk6_ST;
  foo->RKstepper5 = &TSIL_rk5_ST;

  /* Set values in the data object */
  /* DGR -- leave y,z undefined, or set to nan or something? */
  foo->x  = x;
  foo->y  = TSIL_Infinity;
  foo->z  = TSIL_Infinity;
  foo->u  = u;
  foo->v  = v;
  foo->qq = qq;

  /* Set default scale factor */
  foo->scaleFac = 1.0L;

  /* Calculate thresholds... */
/*   foo->threshold[0] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(z), 2); */
/*   foo->threshold[1] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(u), 2); */
  foo->threshold[0] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(u) + TSIL_SQRT(v), 2);
/*   foo->threshold[3] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(z) + TSIL_SQRT(v), 2); */
  foo->threshMin    = foo->threshold[0];
  foo->nThresh = 1;

  /* ...and pseudo-thresholds: */
/*   foo->pseudoThreshold[0] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(z), 2); */
/*   foo->pseudoThreshold[1] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(u), 2); */
  foo->pseudoThreshold[0] = TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(u) - TSIL_SQRT(v), 2);
  foo->pseudoThreshold[1] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(u) + TSIL_SQRT(v), 2);
  foo->pseudoThreshold[2] = TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(u) - TSIL_SQRT(v), 2);
/*   foo->pseudoThreshold[5] = TSIL_POW(TSIL_SQRT(y) + TSIL_SQRT(z) - TSIL_SQRT(v), 2); */
/*   foo->pseudoThreshold[6] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(z) + TSIL_SQRT(v), 2); */
/*   foo->pseudoThreshold[7] = TSIL_POW(TSIL_SQRT(y) - TSIL_SQRT(z) - TSIL_SQRT(v), 2); */
  foo->psThreshMin        = TSIL_MinAbs (foo->pseudoThreshold, 3);
  foo->nPthresh = 3;

  /* Set values in all sub-objects: */
/*   TSIL_ConstructB (&(foo->B[xz]), x, z, qq); */
/*   TSIL_ConstructB (&(foo->B[yu]), y, u, qq); */

/*   TSIL_ConstructS (&(foo->S[vyz]), vyz, v, y, z, qq); */
  TSIL_ConstructS (&(foo->S[uxv]), uxv, u, x, v, qq);

/*   TSIL_ConstructT (&(foo->T[vyz]), vyz, v, y, z, qq); */
  TSIL_ConstructT (&(foo->T[uxv]), uxv, u, x, v, qq);
/*   TSIL_ConstructT (&(foo->T[yzv]), yzv, y, z, v, qq); */
  TSIL_ConstructT (&(foo->T[xuv]), xuv, x, u, v, qq);
/*   TSIL_ConstructT (&(foo->T[zyv]), zyv, z, y, v, qq); */
  TSIL_ConstructT (&(foo->T[vxu]), vxu, v, x, u, qq);

/*   TSIL_ConstructU (&(foo->U[zxyv]), zxyv, z, x, y, v, qq); */
/*   TSIL_ConstructU (&(foo->U[uyxv]), uyxv, u, y, x, v, qq); */
/*   TSIL_ConstructU (&(foo->U[xzuv]), xzuv, x, z, u, v, qq); */
/*   TSIL_ConstructU (&(foo->U[yuzv]), yuzv, y, u, z, v, qq); */

/*   TSIL_ConstructV (&(foo->V[zxyv]), zxyv, z, x, y, v, qq); */
/*   TSIL_ConstructV (&(foo->V[uyxv]), uyxv, u, y, x, v, qq); */
/*   TSIL_ConstructV (&(foo->V[xzuv]), xzuv, x, z, u, v, qq); */
/*   TSIL_ConstructV (&(foo->V[yuzv]), yuzv, y, u, z, v, qq); */

/*   mtype = TSIL_ConstructM (&(foo->M), x, y, z, u, v, qq); */

  /* Set Runge-Kutta step-size parameter parameters.  */
  /* These can be reset by calling TSIL_ResetStepSizeParams */

  foo->precisionGoal = TSIL_PRECISION_GOAL;
  foo->nStepsStart   = TSIL_NSTEPS_START;
  foo->nStepsMin     = TSIL_NSTEPS_MIN;
  foo->nStepsMaxCon  = TSIL_NSTEPS_MAX_CON;
  foo->nStepsMaxVar  = TSIL_NSTEPS_MAX_VAR;

  foo->isInitialized = TRUE;
  foo->status        = UNEVALUATED;

  return 0;
}
