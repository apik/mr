/* ******************************************************************* */
/*              EVALUATION DRIVER PROGRAM                              */
/* ******************************************************************* */

#include "internal.h"
#include "tsil_params.h"

/* ******************************************************************* */

int TSIL_Evaluate (TSIL_DATA *foo, TSIL_REAL s)
{
  int status = 1;
  int isAnalytic = 0;
  TSIL_REAL x0, y0, z0, u0, v0, s0, qq0;
  int tmpWarns;

  if (foo->isInitialized != YES)
    TSIL_Error("TSIL_Evaluate",
	       "You must first set parameter values using TSIL_SetParams!", 3);

  /* Set s value in foo */
  foo->s = s;

  /* Decide if this is a special case known analytically, or a generic
     case that requires Runge-Kutta evaluation: */

  /* Temporarily disable WARNs */
  tmpWarns = printWarns;
  printWarns = NO;

  if (foo->whichFns == STUM)
    isAnalytic = TSIL_Manalytic (foo->x, foo->y, foo->z, foo->u, foo->v, 
				 foo->s, &(foo->M.value));
  else if (foo->whichFns == STU)
    isAnalytic = TSIL_Uanalytic (foo->x, foo->z, foo->u, foo->v, 
				 foo->s, foo->qq, &(foo->U[xzuv].value));
  else if (foo->whichFns == ST) {
    isAnalytic = TSIL_Tanalytic (foo->x, foo->u, foo->v, 
				 foo->s, foo->qq, &(foo->T[xuv].value));
    isAnalytic *= TSIL_Tanalytic (foo->u, foo->x, foo->v, 
				  foo->s, foo->qq, &(foo->T[uxv].value));
    isAnalytic *= TSIL_Tanalytic (foo->v, foo->x, foo->u, 
				  foo->s, foo->qq, &(foo->T[vxu].value));
  }
  else
    TSIL_Error ("TSIL_Evaluate",
		"This can't happen: whichFns was not properly set!", 2);

  /* Restore warnings */
/*   printWarns = YES; */

  if (isAnalytic == TRUE) {
    if (0 == TSIL_CaseSpecial (foo)) 
      TSIL_Error("TSIL_Evaluate",
		 "This can't happen! TSIL_CaseSpecial returned 0 when isAnalytic is TRUE.", 1);
    /* ...and we are finished in this case */
    foo->status = ANALYTIC;
  }
  else {
    /* Generic evaluation by integration: */
    /* In v1.2 we store the initial args and restore them exactly: */
    x0 = foo->x; y0 = foo->y; z0 = foo->z; u0 = foo->u; 
    v0 = foo->v; s0 = s; qq0 = foo->qq;

    TSIL_Rescale (foo);

    /* Check if this is the "unnatural" threshold case */
    if (TSIL_UnnaturalCase (foo))
      TSIL_Warn ("TSIL_Evaluate",
		 "'Unnatural' threshold case! Expect reduced accuracy.");

    /* Generic evaluation */
    TSIL_CaseGeneric (foo);

    if (1 == REDOANALYTIC) {
      /*
	The following line puts in the analytic results for B, and if
	available, S,T,U. But leaving it out allows debugging of the
	RK evaluation of analytic B,S,T,U cases, which is quite useful
	for now.
      */
      TSIL_CaseSpecial (foo);
    }

    /* Undo rescaling */
    TSIL_Unscale (foo);
    foo->x = x0; foo->y = y0; foo->z = z0; foo->u = u0; 
    foo->v = v0; s = s0; foo->qq = qq0;
  }

  /* Set additional functions: */
  TSIL_SetTbar (foo);
  TSIL_SetV (foo);
  TSIL_SetBold (foo);

/*   if (foo->whichFns == STU) */
/*     foo->B[xz].value = B(foo->x, foo->z, foo->s, foo->qq); */

  /* Restore warnings */
  printWarns = tmpWarns;

  /* Implement status codes eventually */
  return status;
}
