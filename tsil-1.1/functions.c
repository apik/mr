/* Miscellaneous and general functions */

#include "internal.h"
#include <string.h>

/* ******************************************************************* */
/* Implements Eq. 2.29 of hep-ph/0307101                               */

TSIL_COMPLEX Delta (TSIL_COMPLEX a, TSIL_COMPLEX b, TSIL_COMPLEX c)
{
  return a*a + b*b + c*c - 2.0L*a*b - 2.0L*a*c - 2.0L*b*c;
}

/* ******************************************************************* */

TSIL_REAL MaxAbs (TSIL_REAL *z, int n)
{
  TSIL_REAL test, res = TSIL_FABS(*z);
  while (n > 1) {
    test = TSIL_FABS(*(z+n-1));
    if (test > res)
      res = test;
    n--;
  }
  return res;
}

/* ******************************************************************* */

TSIL_REAL MinAbs (TSIL_REAL *z, int n)
{
  TSIL_REAL test, res = TSIL_FABS(*z);
  while (n > 1) {
    test = TSIL_FABS(*(z+n-1));
    if (test < res)
      res = test;
    n--;
  }
  return res;
}

/* ******************************************************************* */

TSIL_COMPLEX AddIeps (TSIL_COMPLEX S)
{
  if (TSIL_FABS(TSIL_CIMAG(S))/TSIL_CABS(S) < TSIL_TOL)
    return TSIL_CREAL(S) + TSIL_CABS(S)*I*TSIL_TOL*TSIL_TOL;
  else
    return S;
}

/* ******************************************************************* */

TSIL_COMPLEX EtaBranch (TSIL_COMPLEX z1, TSIL_COMPLEX z2)
{
  return TSIL_CLOG(z1*z2) - TSIL_CLOG(z1) - TSIL_CLOG(z2);
}

/* ******************************************************************* */
/*
  Returns 1 (TRUE) if the evaluation point is within THRESH_CUTOFF of
  a threshold, and sets *sthresh to the nearby threshold value.
  Returns 0 (FALSE) if we are further away than THRESH_CUTOFF from all
  thresholds.
*/

int NearThreshold (TSIL_DATA *foo, TSIL_REAL *sthresh, TSIL_REAL mindistance)
{
  TSIL_REAL distance;
  int i, nThresh;
  int areWeCloseToAThreshold = NO;

  *sthresh = 0.0L; /* This should be ignored if it survives. */

  for (i=0; i<(foo->nThresh); i++) {
    distance = TSIL_FABS((foo->s) - (foo->threshold[i]));
    if (distance < mindistance) {
      mindistance = distance;
      *sthresh = foo->threshold[i];
      areWeCloseToAThreshold = YES;
    }
  }
  return areWeCloseToAThreshold;
}

/* ******************************************************************* */

void ScaleData (TSIL_DATA *foo, TSIL_REAL sfac)
{
  TSIL_REAL sfToThe[5];
  TSIL_REAL extrascalefac1, extrascalefac2;
  int i, j;

  /* Update overall scale factor */
  foo->scaleFac *= sfac;

  sfToThe[0] = 1.0L;
  for (i=1; i<5; i++)
    sfToThe[i] = sfac * sfToThe[i-1];

  /* Now rescale everything as needed... */
  foo->x  /= sfac;
  foo->y  /= sfac;
  foo->z  /= sfac;
  foo->u  /= sfac;
  foo->v  /= sfac;
  foo->s  /= sfac;
  foo->qq /= sfac;

  for (i=0; i<foo->nThresh; i++)
    foo->threshold[i] /= sfac;

  for (i=0; i<foo->nPthresh; i++)
    foo->pseudoThreshold[i] /= sfac;

  foo->threshMin /= sfac;
  foo->psThreshMin /= sfac;

  /* B-type functions */
  for (i=0; i<2; i++) {
    foo->B[i].deriv *= sfac;
    for (j=0; j<2; j++) {
      foo->B[i].arg[j] /= sfac;
      foo->B[i].B_den[j] /= sfToThe[1];
      foo->B[i].B_cB[j] /= sfToThe[1];
      foo->B[i].B_c[j] /= sfToThe[1];
    }
  }

  /* S-type functions */
  for (i=0; i<2; i++) { 
    foo->S[i].value /= sfac;
    foo->S[i].S_c /= sfac;
    for (j=0; j<3; j++)
      foo->S[i].arg[j] /= sfac;
  }

  /* T-type functions */
  for (i=0; i<6; i++) {
    foo->T[i].deriv *= sfac;

    for (j=0; j<3; j++) {
      foo->T[i].arg[j] /= sfac;
      foo->Tbar[i].arg[j] /= sfac;
    }
    for (j=0; j<4; j++) {
      foo->T[i].T_den[j]    /= sfac;
      foo->T[i].cTT1_num[j] /= sfac;
      foo->T[i].cTT2_num[j] /= sfac;
      foo->T[i].cTT3_num[j] /= sfac;
      foo->T[i].cT_num[j] /= sfac;
    }
  }

  /* U-type functions */
  for (i=0; i<4; i++) {
    foo->U[i].deriv *= sfac;

    for (j=0; j<4; j++)
      foo->U[i].arg[j] /= sfac;

    foo->U[i].den_th /= sfToThe[1];
    foo->U[i].den_ps /= sfToThe[1];
    foo->U[i].cU_th /= sfToThe[1]; 
    foo->U[i].cU_ps /= sfToThe[1]; 
    foo->U[i].cT1_th /= sfToThe[1];  
    foo->U[i].cT1_ps /= sfToThe[1]; 
    foo->U[i].cT2 /= sfToThe[1];  
    foo->U[i].cT3 /= sfToThe[1];
    foo->U[i].con /= sfToThe[1];
  }

  /* M-type function */
  for (i=0; i<5; i++)
    foo->M.arg[i] /= sfac;

  foo->M.value *= sfac;

  /* Recall that M.deriv is actually d(sM)/ds, not dM/ds */
  foo->M.deriv *= sfac;

  foo->M.THxz /= sfac;
  foo->M.THyu /= sfac;
  foo->M.PSxz /= sfac;
  foo->M.PSyu /= sfac;

  extrascalefac1 = TSIL_POW(TSIL_SQRT(sfac), foo->M.extramassdim1);
  extrascalefac2 = TSIL_POW(TSIL_SQRT(sfac), foo->M.extramassdim2);

  for (i=0; i<3; i++)
    foo->M.adenom[i] /= sfToThe[3-i]*extrascalefac1;

  for (i=0; i<4; i++) {
    foo->M.bMS[i] /= sfToThe[2]*extrascalefac1;
    foo->M.bM[i] /= sfToThe[3]*extrascalefac1;
    for (j=0; j<2; j++)
      foo->M.bMU[i][j] /= sfToThe[2]*extrascalefac2;
  }

  for (i=0; i<5; i++) {
    foo->M.aMT[i] /= sfToThe[2]*extrascalefac1;
    for (j=0; j<4; j++)
      foo->M.bMT[i][j] /= sfToThe[3]*extrascalefac1;
  }

  for (i=0; i<2; i++) {
    foo->M.aMB[i] /= sfToThe[1]*extrascalefac2;
    for (j=0; j<2; j++)
      foo->M.bMB[i][j] /= sfToThe[2]*extrascalefac2;
  }

  foo->M.aM /= sfToThe[2]*extrascalefac1;
  foo->M.aMT5s /= sfToThe[1];
  foo->M.aMS /= sfToThe[1];
  foo->M.cMSconst /= sfToThe[1];
  foo->M.aMs /= sfToThe[1];

  return;
}

/* ******************************************************************* */
/* Scales all momenta by the largest of x,y,z,u,v,s                    */

void Rescale (TSIL_DATA *foo)
{
  TSIL_REAL tmp[6], sf;
  int i, n = 0;

  /* This is the original code, reorganized to handle different
     cases: */
  tmp[0] = foo->x;
  tmp[1] = foo->u;
  tmp[2] = foo->v;
  tmp[3] = foo->s;
  tmp[4] = foo->z;
  tmp[5] = foo->y;

  if      (foo->whichFns == STUM) n = 6;
  else if (foo->whichFns == STU ) n = 5;
  else if (foo->whichFns == ST  ) n = 4;

  sf = MaxAbs(tmp, n);
  /* End original code */


  /* DGR - pick a more central value? Doesn't seem to work as well as
     using the max value.  Check generation of test data? */

  /* DGR - Reorganized to facilitate handling of different cases below */
/*   if (foo->x > TSIL_TOL) { */
/*     tmp[n] = foo->x; */
/*     n++; */
/*   } */
/*   if (foo->u > TSIL_TOL) { */
/*     tmp[n] = foo->u; */
/*     n++; */
/*   } */
/*   if (foo->v > TSIL_TOL) { */
/*     tmp[n] = foo->v; */
/*     n++; */
/*   } */
/*   if (TSIL_FABS(foo->s) > TSIL_TOL) { */
/*     tmp[n] = foo->s; */
/*     n++; */
/*   } */
/*   if (foo->whichFns != ST && foo->z > TSIL_TOL) { */
/*     tmp[n] = foo->z; */
/*     n++; */
/*   } */
/*   if (foo->whichFns == STUM && foo->y > TSIL_TOL) { */
/*     tmp[n] = foo->y; */
/*     n++; */
/*   } */

  /* Determine the scaling factor */
/*   sf = TSIL_EXP(0.5L*(TSIL_LOG(MaxAbs(tmp, n)) + TSIL_LOG(MinAbs(tmp, n)))); */

  ScaleData (foo, sf);

  return;
}

/* ******************************************************************* */
/* Undoes any rescaling                                                */

void Unscale (TSIL_DATA *foo)
{
  TSIL_REAL sf;

  /* Determine the "unscaling" factor */
  sf = 1.0L/(foo->scaleFac);

  ScaleData (foo, sf);
  foo->scaleFac = 1.0L;

  return;
}

/* ******************************************************************* */

void TSIL_PrintInfo (void)
{
  printf("(* Variable size: ");
#if defined(TSIL_SIZE_DOUBLE)
  printf("DOUBLE");
#else
  printf("LONG DOUBLE");
#endif
  printf(". *)\n");

  return;
}

/* ******************************************************************* */
/* Extract function values from the data object to an array            */

void TSIL_GetData (TSIL_DATA    *foo, 
		   const char   *which, 
		   TSIL_COMPLEX *val)
{
  int i;

  if (foo->status == UNEVALUATED)
    TSIL_Error ("TSIL_GetData", "This case has not been evaluated!", 1);

  if (foo->whichFns == STUM) {
    if (!strcmp(which, "M")) {
      *val = foo->M.value;
    }
    else if (!strcmp(which, "U")) {
      for (i=0; i<NUM_U_FUNCS; i++)
	val[i] = foo->U[i].value;
    }
    else if (!strcmp(which, "T")) {
      for (i=0; i<NUM_T_FUNCS; i++)
	val[i] = foo->T[i].value;
    }
    else if (!strcmp(which, "S")) {
      for (i=0; i<NUM_S_FUNCS; i++)
	val[i] = foo->S[i].value;
    }
    else if (!strcmp(which, "B")) {
      for (i=0; i<NUM_B_FUNCS; i++)
	val[i] = foo->B[i].value;
    }
    else if (!strcmp(which, "V")) {
      for (i=0; i<NUM_V_FUNCS; i++)
	val[i] = foo->V[i].value;
    }
    else if (!strcmp(which, "TBAR")) {
      for (i=0; i<NUM_T_FUNCS; i++)
	val[i] = foo->Tbar[i].value;
    }
    else
      TSIL_Error ("TSIL_GetData", "Invalid identifier", 3);
  }
  else
    TSIL_Error ("TSIL_GetData", 
		"Cannot use this function with subset cases .", 4);

  return;
}

/* ******************************************************************* */
/* Extract "bold" function values from the data object to an array     */

void TSIL_GetBoldData (TSIL_DATA    *foo,
		       const char   *which,
		       TSIL_COMPLEX val[][3])
{
  int i, j;

  if (foo->status == UNEVALUATED)
    TSIL_Error ("TSIL_GetBoldData", "This case has not been evaluated!", 1);

  if (foo->whichFns == STUM) {
    if (!strcmp(which, "U")) {
      for (i=0; i<NUM_U_FUNCS; i++)
	for (j=0; j<3; j++)
	  val[i][j] = foo->U[i].bold[j];
    }
    else if (!strcmp(which, "V")) {
      for (i=0; i<NUM_V_FUNCS; i++)
	for (j=0; j<3; j++)
	  val[i][j] = foo->V[i].bold[j];
    }
    else if (!strcmp(which, "T")) {
      for (i=0; i<NUM_T_FUNCS; i++)
	for (j=0; j<3; j++)
	  val[i][j] = foo->T[i].bold[j];
    }
    else if (!strcmp(which, "S")) {
      for (i=0; i<NUM_S_FUNCS; i++)
	for (j=0; j<3; j++)
	  val[i][j] = foo->S[i].bold[j];
    }
    else
      TSIL_Error ("TSIL_GetBoldData", "Invalid identifier", 3);
  }
  else
    TSIL_Error ("TSIL_GetBoldData",
		"Can't use this for subset evaluation (STU or ST).", 4);

  return;
}

/* ******************************************************************* */
/* Extract a single function value from the data object                */

TSIL_COMPLEX TSIL_GetFunction (TSIL_DATA *foo, const char *which)
{
  /* Below is cut and pasted directly from tsil_names.h */
  const char *uname[] = {"Uzxyv","Uuyxv","Uxzuv","Uyuzv"};
  const char *tname[] = {"Tvyz", "Tuxv", "Tyzv", "Txuv", "Tzyv", "Tvxu"};
  const char *sname[] = {"Svyz", "Suxv"};
  const char *bname[] = {"Bxz", "Byu"};
  const char *vname[] = {"Vzxyv","Vuyxv","Vxzuv","Vyuzv"};
  const char *tbarname[] = {"TBARvyz", "TBARuxv", "TBARyzv", 
                            "TBARxuv", "TBARzyv", "TBARvxu"};
  int i;

  /* Check evaluation status: */
  if (foo->status == UNEVALUATED)
    TSIL_Error ("TSIL_GetFunction",
		"This case has not yet been evaluated!", 1);

  /* For simplicity, do evaluation cases separately */
  if (foo->whichFns == STUM) {
    if (!strncmp(which, "M", 1)) {
      return foo->M.value;
    }
    else if (!strncmp(which, "U", 1)) {
      for (i=0; i<NUM_U_FUNCS; i++)
	if (!strcmp(which, uname[i]))
	  return foo->U[i].value;
    }
    else if (!strncmp(which, "TBAR", 4)) {
      for (i=0; i<NUM_T_FUNCS; i++)
	if (!strcmp(which, tbarname[i]))
	  return foo->Tbar[i].value;
    }
    else if (!strncmp(which, "T", 1)) {
      for (i=0; i<NUM_T_FUNCS; i++)
	if (!strcmp(which, tname[i]))
	  return foo->T[i].value;
    }
    else if (!strncmp(which, "S", 1)) {
      for (i=0; i<NUM_S_FUNCS; i++)
	if (!strcmp(which, sname[i]))
	  return foo->S[i].value;
    }
    else if (!strncmp(which, "B", 1)) {
      for (i=0; i<NUM_B_FUNCS; i++)
	if (!strcmp(which, bname[i]))
	  return foo->B[i].value;
    }
    else if (!strncmp(which, "V", 1)) {
      for (i=0; i<NUM_V_FUNCS; i++)
	if (!strcmp(which, vname[i]))
	  return foo->V[i].value;
    }
    else
      TSIL_Error ("TSIL_GetFunction", "Invalid identifier", 3);
  }
  else if (foo->whichFns == STU) {

    if (!strncmp(which, "U", 1))
	  return foo->U[xzuv].value;

    else if (!strncmp(which, "TBAR", 4)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tbarname[i]))
	  return foo->Tbar[i].value;
    }
    else if (!strncmp(which, "T", 1)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tname[i]))
	  return foo->T[i].value;
    }
    else if (!strncmp(which, "S", 1))
      return foo->S[uxv].value;

    else if (!strncmp(which, "B", 1)) {
      return foo->B[xz].value;
    }
    else if (!strncmp(which, "V", 1))
      return foo->V[xzuv].value;

    else
      TSIL_Error ("TSIL_GetFunction", "Invalid identifier for this case", 3);
  }
  else if (foo->whichFns == ST) {

    if (!strncmp(which, "TBAR", 4)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tbarname[i]))
	  return foo->Tbar[i].value;
    }
    else if (!strncmp(which, "T", 1)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tname[i]))
	  return foo->T[i].value;
    }
    else if (!strncmp(which, "S", 1))
      return foo->S[uxv].value;

/*     else if (!strncmp(which, "B", 1)) { */
/*       for (i=0; i<NUM_B_FUNCS; i++) */
/* 	if (!strcmp(which, bname[i])) */
/* 	  return foo->B[i].value; */
/*     } */
    else
      TSIL_Error ("TSIL_GetFunction", "Invalid identifier for this case", 3);
  }
  else
    TSIL_Error ("TSIL_GetFunction", 
		"This can't happen!  Evaluation mode set incorrectly.", 4);
}

/* ******************************************************************* */
/* Extract a single "bold" function value from the data object         */

TSIL_COMPLEX TSIL_GetBoldFunction (TSIL_DATA  *foo,
				   const char *which,
				   int        n)
{
  /* Below is cut and pasted directly from tsil_names.h */
  const char *uname[] = {"Uzxyv","Uuyxv","Uxzuv","Uyuzv"};
  const char *tname[] = {"Tvyz", "Tuxv", "Tyzv", "Txuv", "Tzyv", "Tvxu"};
  const char *sname[] = {"Svyz", "Suxv"};
  const char *vname[] = {"Vzxyv","Vuyxv","Vxzuv","Vyuzv"};
  int i;

  /* Check inputs and evaluation status: */
  if (n<0 || n>2)
    TSIL_Error ("TSIL_GetBoldFunction",
		"Invalid power specified, must be 0,1, or 2.", 2);

  if (foo->status == UNEVALUATED)
    TSIL_Error ("TSIL_GetBoldFunction",
		"This case has not been evaluated!", 1);

  /* Again, branch on evaluation case */

  if (foo->whichFns == STUM) {
    /* Find and return appropriate value: */
    if (!strncmp(which, "U", 1)) {
      for (i=0; i<NUM_U_FUNCS; i++)
	if (!strcmp(which, uname[i]))
	  return foo->U[i].bold[n];
    }
    else if (!strncmp(which, "T", 1)) {
      for (i=0; i<NUM_T_FUNCS; i++)
	if (!strcmp(which, tname[i]))
	  return foo->T[i].bold[n];
    }
    else if (!strncmp(which, "S", 1)) {
      for (i=0; i<NUM_S_FUNCS; i++)
	if (!strcmp(which, sname[i]))
	  return foo->S[i].bold[n];
    }
    else if (!strncmp(which, "V", 1)) {
      for (i=0; i<NUM_V_FUNCS; i++)
	if (!strcmp(which, vname[i]))
	  return foo->V[i].bold[n];
    }
    else
      TSIL_Error ("TSIL_GetBoldFunction",
		  "Invalid identifier for this case", 3);
  }
  else if (foo->whichFns == STU) {
    
    /* Find and return appropriate value: */
    if (!strncmp(which, "U", 1))
      return foo->U[xzuv].bold[n];

    else if (!strncmp(which, "T", 1)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tname[i]))
	  return foo->T[i].bold[n];
    }
    else if (!strncmp(which, "S", 1))
      return foo->S[uxv].bold[n];

    else if (!strncmp(which, "V", 1))
      return foo->V[xzuv].bold[n];
    
    else
      TSIL_Error ("TSIL_GetBoldFunction", 
		  "Invalid identifier for this case.", 3);
  }
  else if (foo->whichFns == ST) {

    if (!strncmp(which, "T", 1)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tname[i]))
	  return foo->T[i].bold[n];
    }
    else if (!strncmp(which, "S", 1))
      return foo->S[uxv].bold[n];

    else
      TSIL_Error ("TSIL_GetBoldFunction",
		  "Invalid identifier for this case.", 3);
  }
  else
    TSIL_Error ("TSIL_GetBoldFunction", 
		"This can't happen!  Evaluation mode set incorrectly.", 4);
}

/* ******************************************************************* */
/* Returns true/false indicating if z is TSIL_Infinity               */

int TSIL_IsInfinite (TSIL_COMPLEX z)
{
  /* DGR - modified to work with gcc4 */
#if __GNUC__==4
  if (isnan(TSIL_CREAL(z)) || isinf(TSIL_CREAL(z)) ||
      isnan(TSIL_CIMAG(z)) || isinf(TSIL_CIMAG(z)))
    return TRUE;
  else
    return FALSE;
#else
  if (isnan(TSIL_CREAL(z)) || isnan(TSIL_CIMAG(z)))
    return TRUE;
  else
    return FALSE;
#endif
}

/* ******************************************************************* */
/* Generic printing of complexes                                       */

void TSIL_cprintf (TSIL_COMPLEX z)
{
  cfprintf (stdout, (double complex) z);
  return;
}

/* ******************************************************************* */
/* Mathematica-compatible printing of complexes                        */

void TSIL_cprintfM (TSIL_COMPLEX z)
{
  cfprintfM (stdout, (double complex) z);
  return;
}

/* ******************************************************************* */
/* Generic printing of complexes                                       */

void cfprintf (FILE *fp, double complex z)
{
  if (TSIL_IsInfinite (z))
    fprintf(fp, " ComplexInfinity");
  else
    fprintf(fp, "% 18.16lf, % 18.16lf", creal(z), cimag(z));

  return;
}

/* ******************************************************************* */
/* Mathematica-compatible printing of complexes                        */

void cfprintfM (FILE *fp, double complex z)
{
  if (TSIL_IsInfinite (z))
    fprintf(fp, " ComplexInfinity");
  else
    fprintf(fp, "% 18.16lf + % 18.16lf I", creal(z), cimag(z));

  return;
}

/* ******************************************************************* */
/* Returns the intrinsic data size used in building the library        */

int TSIL_DataSize (void)
{
#if defined(TSIL_SIZE_DOUBLE)
  return DOUBLE;
#else
  return LONG_DOUBLE;
#endif
}

/* ******************************************************************* */

void TSIL_PrintData (TSIL_DATA *foo)
{
  TSIL_WriteData (stdout, foo);
  return;
}

/* ******************************************************************* */

void TSIL_PrintDataM (TSIL_DATA *foo)
{
  TSIL_WriteDataM (stdout, foo);
  return;
}

/* ******************************************************************* */

void TSIL_WriteData (FILE *fp, TSIL_DATA *foo)
{
#include "tsil_names.h"

  TSIL_REAL fac = foo->scaleFac;
  int j, k;

  if (foo->status == UNEVALUATED) {
    TSIL_Warn ("Write/TSIL_PrintData", "This case has not yet been evaluated!");
    return;
  }

  /* Format strings shortened somewhat */
  fprintf(fp, "x = %.12Lf\n", (long double) (fac * foo->x));
  if (foo->whichFns == STUM)
    fprintf(fp, "y = %.12Lf\n", (long double) (fac * foo->y));
  if (foo->whichFns != ST)
    fprintf(fp, "z = %.12Lf\n", (long double) (fac * foo->z));
  fprintf(fp, "u = %.12Lf\n", (long double) (fac * foo->u));
  fprintf(fp, "v = %.12Lf\n", (long double) (fac * foo->v));
  fprintf(fp, "s = %.12Lf\n", (long double) (fac * foo->s));
  fprintf(fp, "qq = %.12Lf\n", (long double) (fac * foo->qq));

  fprintf(fp, "\n");

  /* Again branch on evaluation mode for simplicity */
  if (foo->whichFns == STUM) {

    fprintf(fp, "Mxyzuv   = ");
    TSIL_cprintf(foo->M.value); fprintf(fp, "\n");

    fprintf(fp, "\n");
    for (j=0; j<4; j++) {
      fprintf(fp, "%s    = ",uname[j]);
      TSIL_cprintf(foo->U[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<6; j++) {
      fprintf(fp, "%s     = ",tname[j]);
      TSIL_cprintf(foo->T[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<2; j++) {
      fprintf(fp, "%s     = ",sname[j]);
      TSIL_cprintf(foo->S[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<2; j++) {
      fprintf(fp, "%s      = ",bname[j]);
      TSIL_cprintf(foo->B[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<4; j++) {
      fprintf(fp, "%s    = ",vname[j]);
      TSIL_cprintf(foo->V[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<6; j++) {
      fprintf(fp, "%s  = ",tbarname[j]);
      TSIL_cprintf(foo->Tbar[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<4; j++) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s  = ",uuname[j][k]);
	TSIL_cprintf(foo->U[j].bold[k]); fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
    for (j=0; j<4; j++) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s  = ",vvname[j][k]);
	TSIL_cprintf(foo->V[j].bold[k]); fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
    for (j=0; j<6; j++) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s   = ",ttname[j][k]);
	TSIL_cprintf(foo->T[j].bold[k]); fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
    for (j=0; j<2; j++) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s   = ",ssname[j][k]);
	TSIL_cprintf(foo->S[j].bold[k]); fprintf(fp, "\n");
      }
    }
  }
  else if (foo->whichFns == STU) {

    fprintf(fp, "\n");
    fprintf(fp, "%s    = ",uname[xzuv]);
    TSIL_cprintf(foo->U[xzuv].value); fprintf(fp, "\n");

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      fprintf(fp, "%s     = ",tname[j]);
      TSIL_cprintf(foo->T[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fprintf(fp, "%s     = ",sname[uxv]);
    TSIL_cprintf(foo->S[uxv].value); fprintf(fp, "\n");

    fprintf(fp, "\n");
    fprintf(fp, "%s      = ",bname[xz]);
    TSIL_cprintf(foo->B[xz].value); fprintf(fp, "\n");

    fprintf(fp, "\n");
    fprintf(fp, "%s    = ",vname[xzuv]);
    TSIL_cprintf(foo->V[xzuv].value); fprintf(fp, "\n");

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      fprintf(fp, "%s  = ",tbarname[j]);
      TSIL_cprintf(foo->Tbar[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (k=0; k<3; k++) {
      fprintf(fp, "%s  = ",uuname[xzuv][k]);
      TSIL_cprintf(foo->U[xzuv].bold[k]); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (k=0; k<3; k++) {
      fprintf(fp, "%s  = ",vvname[xzuv][k]);
      TSIL_cprintf(foo->V[xzuv].bold[k]); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s   = ",ttname[j][k]);
	TSIL_cprintf(foo->T[j].bold[k]); fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
    for (k=0; k<3; k++) {
      fprintf(fp, "%s   = ",ssname[uxv][k]);
      TSIL_cprintf(foo->S[uxv].bold[k]); fprintf(fp, "\n");
    }
  }
  else if (foo->whichFns == ST) {

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      fprintf(fp, "%s     = ",tname[j]);
      TSIL_cprintf(foo->T[j].value); fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fprintf(fp, "%s     = ",sname[uxv]);
    TSIL_cprintf(foo->S[uxv].value); fprintf(fp, "\n");

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      fprintf(fp, "%s  = ",tbarname[j]);
      TSIL_cprintf(foo->Tbar[j].value); fprintf(fp, "\n");
    }

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s   = ",ttname[j][k]);
	TSIL_cprintf(foo->T[j].bold[k]); fprintf(fp, "\n");
      }
    }
    fprintf(fp, "\n");
    for (k=0; k<3; k++) {
      fprintf(fp, "%s   = ",ssname[uxv][k]);
      TSIL_cprintf(foo->S[uxv].bold[k]); fprintf(fp, "\n");
    }
  }
  return;
}

/* ******************************************************************* */
/* Print output in Mathematica-compatible format                       */

void TSIL_WriteDataM (FILE *fp, TSIL_DATA *foo)
{
#include "tsil_names.h"

  TSIL_REAL fac = foo->scaleFac;
  int j, k;

  if (foo->status == UNEVALUATED) {
    TSIL_Warn ("TSIL_PrintDataM", "This case has not been evaluated!");
    return;
  }

  fprintf(fp, "x = %.12Lf;\n", (long double) (fac * foo->x));
  if (foo->whichFns == STUM)
    fprintf(fp, "y = %.12Lf;\n", (long double) (fac * foo->y));
  if (foo->whichFns != ST)
    fprintf(fp, "z = %.12Lf;\n", (long double) (fac * foo->z));
  fprintf(fp, "u = %.12Lf;\n", (long double) (fac * foo->u));
  fprintf(fp, "v = %.12Lf;\n", (long double) (fac * foo->v));
  fprintf(fp, "s = %.12Lf;\n", (long double) (fac * foo->s));
  fprintf(fp, "qq = %.12Lf;\n", (long double) (fac * foo->qq));

  fprintf(fp, "\n");

  /* Branch on evaluation mode */
  if (foo->whichFns == STUM) {
    fprintf(fp, "Mxyzuv   = ");
    TSIL_cprintfM(foo->M.value); fprintf(fp, ";\n");
    
    fprintf(fp, "\n");
    for (j=0; j<4; j++) {
      fprintf(fp, "%s    = ",uname[j]);
      TSIL_cprintfM(foo->U[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<6; j++) {
      fprintf(fp, "%s     = ",tname[j]);
      TSIL_cprintfM(foo->T[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<2; j++) {
      fprintf(fp, "%s     = ",sname[j]);
      TSIL_cprintfM(foo->S[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<2; j++) {
      fprintf(fp, "%s      = ",bname[j]);
      TSIL_cprintfM(foo->B[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<4; j++) {
      fprintf(fp, "%s    = ",vname[j]);
      TSIL_cprintfM(foo->V[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<6; j++) {
      fprintf(fp, "%s  = ",tbarname[j]);
      TSIL_cprintfM(foo->Tbar[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (j=0; j<4; j++) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s  = ",uuname[j][k]);
	TSIL_cprintfM(foo->U[j].bold[k]); fprintf(fp, ";\n");
      }
    }
    fprintf(fp, "\n");
    for (j=0; j<4; j++) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s  = ",vvname[j][k]);
	TSIL_cprintfM(foo->V[j].bold[k]); fprintf(fp, ";\n");
      }
    }
    fprintf(fp, "\n");
    for (j=0; j<6; j++) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s   = ",ttname[j][k]);
	TSIL_cprintfM(foo->T[j].bold[k]); fprintf(fp, ";\n");
      }
    }
    fprintf(fp, "\n");
    for (j=0; j<2; j++) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s   = ",ssname[j][k]);
	TSIL_cprintfM(foo->S[j].bold[k]); fprintf(fp, ";\n");
      }
    }
  }
  else if (foo->whichFns == STU) {

    fprintf(fp, "\n");
    fprintf(fp, "%s    = ",uname[xzuv]);
    TSIL_cprintfM(foo->U[xzuv].value); fprintf(fp, ";\n");

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      fprintf(fp, "%s     = ",tname[j]);
      TSIL_cprintfM(foo->T[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    fprintf(fp, "%s     = ",sname[uxv]);
    TSIL_cprintfM(foo->S[uxv].value); fprintf(fp, ";\n");

    fprintf(fp, "\n");
    fprintf(fp, "%s      = ",bname[xz]);
    TSIL_cprintfM(foo->B[xz].value); fprintf(fp, ";\n");

    fprintf(fp, "\n");
    fprintf(fp, "%s    = ",vname[xzuv]);
    TSIL_cprintfM(foo->V[xzuv].value); fprintf(fp, ";\n");

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      fprintf(fp, "%s  = ",tbarname[j]);
      TSIL_cprintfM(foo->Tbar[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (k=0; k<3; k++) {
      fprintf(fp, "%s  = ",uuname[xzuv][k]);
      TSIL_cprintfM(foo->U[xzuv].bold[k]); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (k=0; k<3; k++) {
      fprintf(fp, "%s  = ",vvname[xzuv][k]);
      TSIL_cprintfM(foo->V[xzuv].bold[k]); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s   = ",ttname[j][k]);
	TSIL_cprintfM(foo->T[j].bold[k]); fprintf(fp, ";\n");
      }
    }
    fprintf(fp, "\n");
    for (k=0; k<3; k++) {
      fprintf(fp, "%s   = ",ssname[uxv][k]);
      TSIL_cprintfM(foo->S[uxv].bold[k]); fprintf(fp, ";\n");
    }
  }
  else if (foo->whichFns == ST) {

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      fprintf(fp, "%s     = ",tname[j]);
      TSIL_cprintfM(foo->T[j].value); fprintf(fp, ";\n");
    }
    fprintf(fp, "\n");
    fprintf(fp, "%s     = ",sname[uxv]);
    TSIL_cprintfM(foo->S[uxv].value); fprintf(fp, ";\n");

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      fprintf(fp, "%s  = ",tbarname[j]);
      TSIL_cprintfM(foo->Tbar[j].value); fprintf(fp, ";\n");
    }

    fprintf(fp, "\n");
    for (j=1; j<6; j+=2) {
      for (k=0; k<3; k++) {
	fprintf(fp, "%s   = ",ttname[j][k]);
	TSIL_cprintfM(foo->T[j].bold[k]); fprintf(fp, ";\n");
      }
    }
    fprintf(fp, "\n");
    for (k=0; k<3; k++) {
      fprintf(fp, "%s   = ",ssname[uxv][k]);
      TSIL_cprintfM(foo->S[uxv].bold[k]); fprintf(fp, ";\n");
    }
  }

  return;
}

/* ******************************************************************* */
/* Two body threshold */

TSIL_REAL Th2 (TSIL_REAL x, TSIL_REAL y)
{
  return TSIL_POW(TSIL_SQRT(x) + TSIL_SQRT(y), 2);
}

/* ******************************************************************* */
/* Two body pseudo-threshold */

TSIL_REAL Ps2 (TSIL_REAL x, TSIL_REAL y)
{
  /* DGR - added */
  if (TSIL_FABS(x - y) < TSIL_TOL)
    return 0.0;
  else
    return TSIL_POW(TSIL_SQRT(x) - TSIL_SQRT(y), 2);
}

/* ******************************************************************* */
/* Alpha \equiv A(x)/Sqrt[x] */

TSIL_REAL Alpha (TSIL_REAL x, TSIL_REAL qq)
{
  if (x > TSIL_TOL)
    return A(x, qq)/TSIL_SQRT(x);
  else
    return 0.0L;
}

/* ******************************************************************* */

int TSIL_GetStatus (TSIL_DATA *foo)
{
  return foo->status;
}

/* ******************************************************************* */
/* Prints the evaluation status of the specified data object: whether
   unevaluated, or evaluaed analytically, numerically by integration
   along the real s axis, or numerically by integration along the
   displaced contour. */

void TSIL_PrintStatus (TSIL_DATA *foo)
{
  if (foo->status == UNEVALUATED)
    TSIL_Warn("TSIL_PrintStatus", "Functions not yet evaluated!");
  else {
    printf("(* Evaluation method: ");
    if (foo->status == ANALYTIC)
      printf("Analytic ");
    else if (foo->status == REAXIS)
      printf("Integration along real s axis ");
    else if (foo->status == CONTOUR)
      printf("Integration along displaced contour ");

    printf("*)\n");
  }
  return;
}

/* ******************************************************************* */
/* Returns true/false indicating whether these parameters correspond   */
/* to the "unnatural" threshold case.                                  */

/* DGR - This should work with the subsidiary cases, as the unused
   parameters in these have been set to TSIL_Infinity */

int UnnaturalCase (TSIL_DATA *foo)
{
  TSIL_REAL x, y, z, u, v;
  TSIL_COMPLEX s;
  TSIL_REAL cutoff = 0.0001L;

  /* For convenience: */
  x = foo->x;
  y = foo->y;
  z = foo->z;
  u = foo->u;
  v = foo->v;
  s = foo->s;

  if (x+TSIL_FABS(y-v)+TSIL_FABS(y-Th2(z,u))+TSIL_CABS(s-z) < cutoff ||
      x+TSIL_FABS(y-v)+TSIL_FABS(y-Ps2(z,u))+TSIL_CABS(s-z) < cutoff ||
      y+TSIL_FABS(x-v)+TSIL_FABS(x-Th2(u,z))+TSIL_CABS(s-u) < cutoff ||
      y+TSIL_FABS(x-v)+TSIL_FABS(x-Ps2(u,z))+TSIL_CABS(s-u) < cutoff ||
      z+TSIL_FABS(u-v)+TSIL_FABS(u-Th2(x,y))+TSIL_CABS(s-x) < cutoff ||
      z+TSIL_FABS(u-v)+TSIL_FABS(u-Ps2(x,y))+TSIL_CABS(s-x) < cutoff ||
      u+TSIL_FABS(z-v)+TSIL_FABS(z-Th2(y,x))+TSIL_CABS(s-y) < cutoff ||
      u+TSIL_FABS(z-v)+TSIL_FABS(z-Ps2(y,x))+TSIL_CABS(s-y) < cutoff)
    return 1;
  else
    return 0;
}

/* ******************************************************************* */

void TSIL_PrintVersion (void)
{
  printf("(* TSIL Version: %s *)\n", TSIL_VERSION);
  return;
}

/* ******************************************************************* */

void TSIL_Error (char *func, char *msg, int val)
{
  fprintf (stderr, "ERROR (%s): %s\n", func, msg);
  fprintf (stderr, "Exiting!\n");
  fflush (stdout);
  exit (val);
}

/* ******************************************************************* */

void TSIL_Warn (char *func, char *msg)
{
  fprintf (stderr, "WARNING (%s): %s\n", func, msg);
  fflush (stdout);
  return;
}

/* ******************************************************************* */
/* Currently disabled! */

void TSIL_Info (char *msg)
{
  /*
  fprintf(stderr, "INFO: ");
  fprintf(stderr, msg);
  fprintf(stderr, "\n");
  fflush(stdout);
  */
  return;
}

/* **************************************************************** */

void CheckConsistent (TSIL_COMPLEX arg1, TSIL_COMPLEX arg2)     
{
  /* If they're both nan, that's OK. */
/*   if ((isnan (TSIL_CREAL (arg1)) || isnan (TSIL_CIMAG (arg1)) ) &&  */
/*       (isnan (TSIL_CREAL (arg2)) || isnan (TSIL_CIMAG (arg2)) ))   */
/*     return; */
  
  /* DGR */
  if (TSIL_IsInfinite (arg1) && TSIL_IsInfinite (arg2))
    return;

  /* If they're both zero, that's OK. */
  if (TSIL_CABS (arg1) + TSIL_CABS (arg2) < TSIL_TOL) return;
  
  /* Otherwise they had better agree exactly. */
  if (TSIL_CABS(1.0L - arg1/arg2) > TSIL_TOL)  
    {
      printf("\n *** Failed consistency check. This can NEVER happen. ***\n");
      exit(1);
    }

  return;
}

/* **************************************************************** */

void TSIL_ResetStepSizeParams (TSIL_DATA *foo,
                               TSIL_REAL precisionGoal,
                               int nstepsStart,
                               int nstepsMaxCon,
                               int nstepsMaxVar,
                               int nstepsMin)
{
  foo->precisionGoal = precisionGoal;
  foo->nStepsStart   = nstepsStart;
  foo->nStepsMaxCon  = nstepsMaxCon;
  foo->nStepsMaxVar  = nstepsMaxVar;
  foo->nStepsMin     = nstepsMin;

  return;
}
