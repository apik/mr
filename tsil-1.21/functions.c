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
  int i;
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
  int n = 0;

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

int TSIL_GetData (TSIL_DATA    *foo, 
		  const char   *which, 
		  TSIL_COMPLEX *val)
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
  int retval = 1;
  char funcname[] = "TSIL_GetData";
  char errmsg0[45] = "Function not defined for these parameters: ";
  char errmsg[55];

  strcpy (errmsg, errmsg0);

  if (foo->status == UNEVALUATED) {
    TSIL_Warn (funcname, "This case has not been evaluated!");
    return 0;
  }

  /* DGR - Only makes sense to extract Bs using this function for
     cases other than STUM */

  if (foo->whichFns == STUM) {
    if (!strcmp(which, "M")) {
      if (TSIL_IsInfinite (*val = foo->M.value))
	TSIL_Warn (funcname, strncat (errmsg, "M", 1));
    }
    else if (!strcmp(which, "U")) {
      for (i=0; i<NUM_U_FUNCS; i++)
	if (TSIL_IsInfinite (val[i] = foo->U[i].value)) {
	  TSIL_Warn (funcname, strncat (errmsg, uname[i], 5));	  
	  strcpy (errmsg, errmsg0);
	}
    }
    else if (!strcmp(which, "T")) {
      for (i=0; i<NUM_T_FUNCS; i++)
	if (TSIL_IsInfinite (val[i] = foo->T[i].value)) {
	  TSIL_Warn (funcname, strncat (errmsg, tname[i], 4));
	  strcpy (errmsg, errmsg0);
	}
    }
    else if (!strcmp(which, "S")) {
      for (i=0; i<NUM_S_FUNCS; i++)
	if (TSIL_IsInfinite (val[i] = foo->S[i].value)) {
	  TSIL_Warn (funcname, strncat (errmsg, sname[i], 4));
	  strcpy (errmsg, errmsg0);
	}
    }
    else if (!strcmp(which, "B")) {
      for (i=0; i<NUM_B_FUNCS; i++)
	if (TSIL_IsInfinite (val[i] = foo->B[i].value)) {
	  TSIL_Warn (funcname, strncat (errmsg, bname[i], 3));
	  strcpy (errmsg, errmsg0);
	}
    }
    else if (!strcmp(which, "V")) {
      for (i=0; i<NUM_V_FUNCS; i++)
	if (TSIL_IsInfinite (val[i] = foo->V[i].value)) {
	  TSIL_Warn (funcname, strncat (errmsg, vname[i], 5));
	  strcpy (errmsg, errmsg0);
	}
    }
    else if (!strcmp(which, "TBAR")) {
      for (i=0; i<NUM_T_FUNCS; i++)
	if (TSIL_IsInfinite (val[i] = foo->Tbar[i].value)) {
	  TSIL_Warn (funcname, strncat (errmsg, tbarname[i], 7));
	  strcpy (errmsg, errmsg0);
	}
    }
    else {
      TSIL_Error (funcname, "Invalid identifier", 1);
      retval = 0;
    }
  }
  else {
    /* Either STU or ST evaluation... */
    if (!strcmp(which, "B")) {
      for (i=0; i<NUM_B_FUNCS; i++)
	if (TSIL_IsInfinite (val[i] = foo->B[i].value)) {
	  TSIL_Warn (funcname, strncat (errmsg, bname[i], 3));
	  strcpy (errmsg, errmsg0);
	}
    }
    else {
      TSIL_Error (funcname,
	"Can only use this to extract B functions in subsidiary cases.", 1);
      retval = 1;
    }
  }

  return retval;
}

/* ******************************************************************* */
/* Extract "bold" function values from the data object to an array     */

int TSIL_GetBoldData (TSIL_DATA    *foo,
		      const char   *which,
		      TSIL_COMPLEX val[][3])
{
  /* Below is cut and pasted directly from tsil_names.h */
  const char *uname[] = {"Uzxyv","Uuyxv","Uxzuv","Uyuzv"};
  const char *tname[] = {"Tvyz", "Tuxv", "Tyzv", "Txuv", "Tzyv", "Tvxu"};
  const char *sname[] = {"Svyz", "Suxv"};
  const char *vname[] = {"Vzxyv","Vuyxv","Vxzuv","Vyuzv"};

  int i, j;
  int retval = 1;
  char funcname[] = "TSIL_GetBoldData";
  char errmsg0[45] = "Function not defined for these parameters: ";
  char errmsg[55];

  if (foo->status == UNEVALUATED) {
    TSIL_Warn (funcname, "This case has not been evaluated!");
    return 0;
  }

  /* Again, only get the whole array for the Bs if STU or ST... */
  if (foo->whichFns == STUM) {
    if (!strcmp(which, "U")) {
      for (i=0; i<NUM_U_FUNCS; i++)
	for (j=0; j<3; j++)
	  if (TSIL_IsInfinite (val[i][j] = foo->U[i].bold[j])) {
	    TSIL_Warn (funcname, strncat (errmsg, uname[i], 5));
	    strcpy (errmsg, errmsg0);
	  }
    }
    else if (!strcmp(which, "V")) {
      for (i=0; i<NUM_V_FUNCS; i++)
	for (j=0; j<3; j++)
	  if (TSIL_IsInfinite (val[i][j] = foo->V[i].bold[j])) {
	    TSIL_Warn (funcname, strncat (errmsg, vname[i], 5));
	    strcpy (errmsg, errmsg0);
	  }
    }
    else if (!strcmp(which, "T")) {
      for (i=0; i<NUM_T_FUNCS; i++)
	for (j=0; j<3; j++)
	  if (TSIL_IsInfinite (val[i][j] = foo->T[i].bold[j])) {
	    TSIL_Warn (funcname, strncat (errmsg, tname[i], 4));
	    strcpy (errmsg, errmsg0);
	  }
    }
    else if (!strcmp(which, "S")) {
      for (i=0; i<NUM_S_FUNCS; i++)
	for (j=0; j<3; j++)
	  if (TSIL_IsInfinite (val[i][j] = foo->S[i].bold[j])) {
	    TSIL_Warn (funcname, strncat (errmsg, sname[i], 4));
	    strcpy (errmsg, errmsg0);
	  }
    }
    else {
      TSIL_Error (funcname, "Invalid identifier", 1);
      retval = 0;
    }
  }
  else {
    TSIL_Error (funcname,
	"Can't use this for subsidiary evaluation (STU or ST).", 1);
    retval = 1;
  }

  return retval;
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
  TSIL_COMPLEX result = (TSIL_COMPLEX) 0.0;
  char funcname[] = "TSIL_GetFunction";
  char errmsg[55] = "Function not defined for these parameters: ";

  /* Check evaluation status: */
  if (foo->status == UNEVALUATED) {
    TSIL_Warn (funcname, "This case has not yet been evaluated!");
    return result;
  }

  /* For simplicity, do evaluation cases separately */
  if (foo->whichFns == STUM) {
    if (!strncmp(which, "M", 1)) {
      result = foo->M.value;
/*       return foo->M.value; */
    }
    else if (!strncmp(which, "U", 1)) {
      for (i=0; i<NUM_U_FUNCS; i++)
	if (!strcmp(which, uname[i]))
	  result = foo->U[i].value;
/* 	  return foo->U[i].value; */
    }
    else if (!strncmp(which, "TBAR", 4)) {
      for (i=0; i<NUM_T_FUNCS; i++)
	if (!strcmp(which, tbarname[i]))
	  result = foo->Tbar[i].value;
/* 	  return foo->Tbar[i].value; */
    }
    else if (!strncmp(which, "T", 1)) {
      for (i=0; i<NUM_T_FUNCS; i++)
	if (!strcmp(which, tname[i]))
	  result = foo->T[i].value;
/* 	  return foo->T[i].value; */
    }
    else if (!strncmp(which, "S", 1)) {
      for (i=0; i<NUM_S_FUNCS; i++)
	if (!strcmp(which, sname[i]))
	  result = foo->S[i].value;
/* 	  return foo->S[i].value; */
    }
    else if (!strncmp(which, "B", 1)) {
      for (i=0; i<NUM_B_FUNCS; i++)
	if (!strcmp(which, bname[i]))
	  result = foo->B[i].value;
/* 	  return foo->B[i].value; */
    }
    else if (!strncmp(which, "V", 1)) {
      for (i=0; i<NUM_V_FUNCS; i++)
	if (!strcmp(which, vname[i]))
	  result = foo->V[i].value;
/* 	  return foo->V[i].value; */
    }
    else {
      /* If we get here the identifier was not recognized: */
      TSIL_Error (funcname, "Invalid identifier", 1);
      return (TSIL_COMPLEX) 0.0;
    }
  }
  else if (foo->whichFns == STU) {

    if (!strcmp(which, uname[xzuv]))
      result = foo->U[xzuv].value;
/* 	  return foo->U[xzuv].value; */

    else if (!strncmp(which, "TBAR", 4)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tbarname[i]))
	  result = foo->Tbar[i].value;
/* 	  return foo->Tbar[i].value; */
    }
    else if (!strncmp(which, "T", 1)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tname[i]))
	  result = foo->T[i].value;
/* 	  return foo->T[i].value; */
    }
    else if (!strcmp(which, sname[uxv]))
      result = foo->S[uxv].value;
/*       return foo->S[uxv].value; */

    else if (!strncmp(which, "B", 1)) {
      for (i=0; i<NUM_B_FUNCS; i++)
	if (!strcmp(which, bname[i]))
	  result = foo->B[i].value;
/* 	  return foo->B[i].value; */
    }
    else if (!strcmp(which, vname[xzuv]))
      result = foo->V[xzuv].value;
/* 	  return foo->V[xzuv].value; */

    else {
      /* If we get here the identifier was not recognized or was wrong: */
      TSIL_Error (funcname, "Invalid identifier for this case", 1);
      result = (TSIL_COMPLEX) 0.0;
    }
  }
  else if (foo->whichFns == ST) {
    
    if (!strncmp(which, "TBAR", 4)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tbarname[i]))
	  result = foo->Tbar[i].value;
/* 	  return foo->Tbar[i].value; */
    }
    else if (!strncmp(which, "T", 1)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tname[i]))
	  result = foo->T[i].value;
/* 	  return foo->T[i].value; */
    }
    else if (!strcmp(which, sname[uxv]))
      result = foo->S[uxv].value;
/*       return foo->S[uxv].value; */

    else if (!strncmp(which, "B", 1)) {
      for (i=0; i<NUM_B_FUNCS; i++)
	if (!strcmp(which, bname[i]))
	  result = foo->B[i].value;
/* 	  return foo->B[i].value; */
    }
    else {
      /* If we get here the identifier was not recognized or was wrong: */
      TSIL_Error ("TSIL_GetFunction", "Invalid identifier for this case", 1);
      result = (TSIL_COMPLEX) 0.0;
    }
  }

  if (TSIL_IsInfinite (result))
    TSIL_Warn (funcname, strncat (errmsg, which, 7));

  return result;
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
  TSIL_COMPLEX result = (TSIL_COMPLEX) 0.0;
  char funcname[] = "TSIL_GetBoldFunction";
  char errmsg[55] = "Function not defined for these parameters: ";

  /* Check inputs and evaluation status: */
  if (n<0 || n>2) {
    TSIL_Error (funcname, "Invalid power specified, must be 0,1, or 2.", 1);
    return result;
  }
  if (foo->status == UNEVALUATED) {
    TSIL_Warn (funcname, "This case has not been evaluated!");
    return result;
  }

  /* Again, branch on evaluation case */

  if (foo->whichFns == STUM) {
    /* Find and return appropriate value: */
    if (!strncmp(which, "U", 1)) {
      for (i=0; i<NUM_U_FUNCS; i++)
	if (!strcmp(which, uname[i]))
	  result = foo->U[i].bold[n];
/* 	  return foo->U[i].bold[n]; */
    }
    else if (!strncmp(which, "T", 1)) {
      for (i=0; i<NUM_T_FUNCS; i++)
	if (!strcmp(which, tname[i]))
	  result = foo->T[i].bold[n];
/* 	  return foo->T[i].bold[n]; */
    }
    else if (!strncmp(which, "S", 1)) {
      for (i=0; i<NUM_S_FUNCS; i++)
	if (!strcmp(which, sname[i]))
	  result = foo->S[i].bold[n];
/* 	  return foo->S[i].bold[n]; */
    }
    else if (!strncmp(which, "V", 1)) {
      for (i=0; i<NUM_V_FUNCS; i++)
	if (!strcmp(which, vname[i]))
	  result = foo->V[i].bold[n];
/* 	  return foo->V[i].bold[n]; */
    }
    else {
      /* If we get here the identifier was not recognized: */
      TSIL_Error (funcname, "Invalid identifier for this case", 1);
      result = (TSIL_COMPLEX) 0.0;
    }
  }
  else if (foo->whichFns == STU) {
    
    /* Find and return appropriate value: */
    if (!strcmp(which, uname[xzuv]))
      result = foo->U[xzuv].bold[n];
/*       return foo->U[xzuv].bold[n]; */

    else if (!strncmp(which, "T", 1)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tname[i]))
	  result = foo->T[i].bold[n];
/* 	  return foo->T[i].bold[n]; */
    }
    else if (!strcmp(which, sname[uxv]))
      result = foo->S[uxv].bold[n];
/*       return foo->S[uxv].bold[n]; */

    else if (!strcmp(which, vname[xzuv]))
      result = foo->V[xzuv].bold[n];
/*       return foo->V[xzuv].bold[n]; */
    
    else {
      /* If we get here the identifier was not recognized: */
      TSIL_Error (funcname, "Invalid identifier for this case.", 1);
      result = (TSIL_COMPLEX) 0.0;
    }
  }
  else if (foo->whichFns == ST) {

    if (!strncmp(which, "T", 1)) {
      for (i=1; i<NUM_T_FUNCS; i+=2)
	if (!strcmp(which, tname[i]))
	  result = foo->T[i].bold[n];
/* 	  return foo->T[i].bold[n]; */
    }
    else if (!strcmp(which, sname[uxv]))
      result = foo->S[uxv].bold[n];
/*       return foo->S[uxv].bold[n]; */

    else {
      /* If we get here the identifier was not recognized: */
      TSIL_Error (funcname, "Invalid identifier for this case.", 1);
      result = (TSIL_COMPLEX) 0.0;
    }
  }

  if (TSIL_IsInfinite (result))
    TSIL_Warn (funcname, strncat (errmsg, which, 7));

  return result;
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
#include <execinfo.h>

void TSIL_Error (char *func, char *msg, int val)
{
  void *callstack[128];
  int i, frames = backtrace(callstack, 128);
  char **strs = backtrace_symbols(callstack, frames);

  fprintf (stderr, "ERROR (%s): %s\n", func, msg);

  for (i=0; i<frames; i++)
    fprintf(stderr, "%s\n", strs[i]);

  fflush (stderr);
  exit (val);
}

/* ******************************************************************* */

void TSIL_Warn (char *func, char *msg)
{
  if (printWarns) {
    fprintf (stderr, "WARNING (%s): %s\n", func, msg);
    fflush (stderr);
  }
  return;
}

/* ******************************************************************* */
/* Currently disabled! */

void TSIL_Info (char *msg)
{
/*   fprintf(stderr, "INFO: ");  */
/*   fprintf(stderr, msg);  */
/*   fprintf(stderr, "\n");  */
/*   fflush(stdout);  */
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

/* **************************************************************** */
/* Copies the results in a TSIL_DATA struct to a TSIL_RESULT struct */

int TSIL_CopyResult (TSIL_DATA *foo, TSIL_RESULT *bar)
{
  if (foo->status == UNEVALUATED) {
    TSIL_Warn ("TSIL_CopyResult", 
	       "TSIL data has not yet been evaluated!");
    return 1;
  }

  if (foo->whichFns == 0) {

    bar->x = foo->x;
    bar->y = foo->y;
    bar->z = foo->z;
    bar->u = foo->u;
    bar->v = foo->v;
    bar->s = foo->s;
    bar->qq = foo->qq;
    
    bar->M = TSIL_GetFunction (foo, "M");
    TSIL_GetData (foo, "U", bar->U);
    TSIL_GetData (foo, "V", bar->V);
    TSIL_GetData (foo, "T", bar->T);
    TSIL_GetData (foo, "S", bar->S);
    TSIL_GetData (foo, "B", bar->B);
    TSIL_GetData (foo, "TBAR", bar->TBAR);
  }
  else if (foo->whichFns == 1) {

    bar->x = foo->x;
/*     bar->y = TSIL_Infinity; */
    bar->z = foo->z;
    bar->u = foo->u;
    bar->v = foo->v;
    bar->s = foo->s;
    bar->qq = foo->qq;

    bar->U[xzuv] = TSIL_GetFunction (foo, "Uxzuv");
    bar->V[xzuv] = TSIL_GetFunction (foo, "Vxzuv");
    bar->T[uxv] = TSIL_GetFunction (foo, "Tuxv");
    bar->T[xuv] = TSIL_GetFunction (foo, "Txuv");
    bar->T[vxu] = TSIL_GetFunction (foo, "Tvxu");
    bar->S[uxv] = TSIL_GetFunction (foo, "Suxv");
    bar->B[xz] = TSIL_GetFunction (foo, "Bxz");
    bar->TBAR[uxv] = TSIL_GetFunction (foo, "TBARuxv");
    bar->TBAR[xuv] = TSIL_GetFunction (foo, "TBARxuv");
    bar->TBAR[vxu] = TSIL_GetFunction (foo, "TBARvxu");
  }
  else if (foo->whichFns == 2) {

    bar->x = foo->x;
/*     bar->y = TSIL_Infinity; */
/*     bar->z = TSIL_Infinity; */
    bar->u = foo->u;
    bar->v = foo->v;
    bar->s = foo->s;
    bar->qq = foo->qq;

    bar->T[uxv] = TSIL_GetFunction (foo, "Tuxv");
    bar->T[xuv] = TSIL_GetFunction (foo, "Txuv");
    bar->T[vxu] = TSIL_GetFunction (foo, "Tvxu");
    bar->S[uxv] = TSIL_GetFunction (foo, "Suxv");
    bar->TBAR[uxv] = TSIL_GetFunction (foo, "TBARuxv");
    bar->TBAR[xuv] = TSIL_GetFunction (foo, "TBARxuv");
    bar->TBAR[vxu] = TSIL_GetFunction (foo, "TBARvxu");
  }
  else 
    TSIL_Error ("TSIL_CopyResult",
		"This can't happen! whichFns was set incorrectly",
		99);

  return 0;
}

/* **************************************************************** */
/* Permutes a TSIL_RESULT.  Allowed permutations are:               */
/*   0) None; just copies  (which = 0 or NOSWAP)                    */
/*   1) x<->y and z<->u    (which = 1 or XYandZU)                   */
/*   2) x<->z and y<->u    (which = 2 or XZandYU)                   */
/*   3) x<->u and y<->z    (which = 3 or XUandYZ)                   */

void TSIL_PermuteResult (TSIL_RESULT *in, 
			 int which, 
			 TSIL_RESULT *out)
{
  /* Always the same: */
  out->v  = in->v;
  out->s  = in->s;
  out->qq = in->qq;
  out->M  = in->M;

  if (which == XYandZU) {

    out->x = in->y;
    out->y = in->x;
    out->z = in->u;
    out->u = in->z;

    out->U[uyxv] = in->U[zxyv];
    out->U[zxyv] = in->U[uyxv];
    out->U[yuzv] = in->U[xzuv];
    out->U[xzuv] = in->U[yuzv];

    out->V[uyxv] = in->V[zxyv];
    out->V[zxyv] = in->V[uyxv];
    out->V[yuzv] = in->V[xzuv];
    out->V[xzuv] = in->V[yuzv];

    out->T[vxu] = in->T[vyz];
    out->T[zyv] = in->T[uxv];
    out->T[xuv] = in->T[yzv];
    out->T[yzv] = in->T[xuv];
    out->T[uxv] = in->T[zyv];
    out->T[vyz] = in->T[vxu];

    out->S[uxv] = in->S[vyz];
    out->S[vyz] = in->S[uxv];

    out->B[yu] = in->B[xz];
    out->B[xz] = in->B[yu];

    out->TBAR[vxu] = in->TBAR[vyz];
    out->TBAR[zyv] = in->TBAR[uxv];
    out->TBAR[xuv] = in->TBAR[yzv];
    out->TBAR[yzv] = in->TBAR[xuv];
    out->TBAR[uxv] = in->TBAR[zyv];
    out->TBAR[vyz] = in->TBAR[vxu];
  }
  else if (which == XZandYU) {

    out->x = in->z;
    out->z = in->x;
    out->y = in->u;
    out->u = in->y;

    out->U[xzuv] = in->U[zxyv];
    out->U[yuzv] = in->U[uyxv];
    out->U[zxyv] = in->U[xzuv];
    out->U[uyxv] = in->U[yuzv];

    out->V[xzuv] = in->V[zxyv];
    out->V[yuzv] = in->V[uyxv];
    out->V[zxyv] = in->V[xzuv];
    out->V[uyxv] = in->V[yuzv];

    out->T[vxu] = in->T[vyz];
    out->T[yzv] = in->T[uxv];
    out->T[uxv] = in->T[yzv];
    out->T[zyv] = in->T[xuv];
    out->T[xuv] = in->T[zyv];
    out->T[vyz] = in->T[vxu];

    out->S[uxv] = in->S[vyz];
    out->S[vyz] = in->S[uxv];

    out->B[yu] = in->B[yu];
    out->B[xz] = in->B[xz];

    out->TBAR[vxu] = in->TBAR[vyz];
    out->TBAR[yzv] = in->TBAR[uxv];
    out->TBAR[uxv] = in->TBAR[yzv];
    out->TBAR[zyv] = in->TBAR[xuv];
    out->TBAR[xuv] = in->TBAR[zyv];
    out->TBAR[vyz] = in->TBAR[vxu];
  }
  else if (which == XUandYZ) {

    out->x = in->u;
    out->u = in->x;
    out->z = in->y;
    out->y = in->z;

    out->U[yuzv] = in->U[zxyv];
    out->U[xzuv] = in->U[uyxv];
    out->U[uyxv] = in->U[xzuv];
    out->U[zxyv] = in->U[yuzv];

    out->V[yuzv] = in->V[zxyv];
    out->V[xzuv] = in->V[uyxv];
    out->V[uyxv] = in->V[xzuv];
    out->V[zxyv] = in->V[yuzv];

    out->T[vyz] = in->T[vyz];
    out->T[xuv] = in->T[uxv];
    out->T[zyv] = in->T[yzv];
    out->T[uxv] = in->T[xuv];
    out->T[yzv] = in->T[zyv];
    out->T[vxu] = in->T[vxu];

    out->S[uxv] = in->S[uxv];
    out->S[vyz] = in->S[vyz];

    out->B[yu] = in->B[xz];
    out->B[xz] = in->B[yu];

    out->TBAR[vyz] = in->TBAR[vyz];
    out->TBAR[xuv] = in->TBAR[uxv];
    out->TBAR[zyv] = in->TBAR[yzv];
    out->TBAR[uxv] = in->TBAR[xuv];
    out->TBAR[yzv] = in->TBAR[zyv];
    out->TBAR[vxu] = in->TBAR[vxu];
  }
  else if (which == NOSWAP) {

    out->x = in->x;
    out->y = in->y;
    out->z = in->z;
    out->u = in->u;

    out->U[uyxv] = in->U[uyxv];
    out->U[zxyv] = in->U[zxyv];
    out->U[yuzv] = in->U[yuzv];
    out->U[xzuv] = in->U[xzuv];

    out->V[uyxv] = in->V[uyxv];
    out->V[zxyv] = in->V[zxyv];
    out->V[yuzv] = in->V[yuzv];
    out->V[xzuv] = in->V[xzuv];

    out->T[vxu] = in->T[vxu];
    out->T[zyv] = in->T[zyv];
    out->T[xuv] = in->T[xuv];
    out->T[yzv] = in->T[yzv];
    out->T[uxv] = in->T[uxv];
    out->T[vyz] = in->T[vyz];

    out->S[uxv] = in->S[uxv];
    out->S[vyz] = in->S[vyz];

    out->B[yu] = in->B[yu];
    out->B[xz] = in->B[xz];

    out->TBAR[vxu] = in->TBAR[vxu];
    out->TBAR[zyv] = in->TBAR[zyv];
    out->TBAR[xuv] = in->TBAR[xuv];
    out->TBAR[yzv] = in->TBAR[yzv];
    out->TBAR[uxv] = in->TBAR[uxv];
    out->TBAR[vyz] = in->TBAR[vyz];
  }
  else
    TSIL_Error ("TSIL_PermuteResults",
		"Invalid permutation specified",
		9);

  return;
}
