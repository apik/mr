/* 
   === fevaluate.c ===

   This is a simple wrapper routine that provides an interface to
   Fortran.  In the Fortran program, include the header
   'tsil_fort.inc'; this defines a common block corresponding to the
   struct below.  One then calls the subroutine

   CALL tsilfevaluate(...)

   defined herein, to evaluate the basis integrals.  

   This code can serve as a model if users are interested in writing
   their own variants with more specialized capabilities.  See the
   file README.txt for detailed notes on the issues involved.
*/

#include "internal.h"
#include "tsil_fortran.h"

/* This struct is used to pass the results of the calculation back to
   the calling program; it matches the data COMMON block defined in
   tsil_fort.inc.  Note that the type sizes MUST agree between them;
   the struct and COMMON block must have the same memory footprint. */

struct {

  /* Basic functions */
  TSIL_COMPLEX_F M;
  TSIL_COMPLEX_F U[NUM_U_FUNCS];
  TSIL_COMPLEX_F T[NUM_T_FUNCS];
  TSIL_COMPLEX_F S[NUM_S_FUNCS];
  TSIL_COMPLEX_F B[NUM_B_FUNCS];
  TSIL_COMPLEX_F V[NUM_V_FUNCS];
  TSIL_COMPLEX_F Tbar[NUM_T_FUNCS];

  /* Bold functions */
  TSIL_COMPLEX_F UU[NUM_U_FUNCS][3];
  TSIL_COMPLEX_F VV[NUM_V_FUNCS][3];
  TSIL_COMPLEX_F TT[NUM_T_FUNCS][3];
  TSIL_COMPLEX_F SS[NUM_S_FUNCS][3];

} results_;

/* For check that data types match (not foolproof!) */
struct {
  int rsize;
} rsize_;

/* ******************************************************************* */
/* Generic wrapper function                                            */

void tsilfevaluate_ (TSIL_REAL_F *x, 
		     TSIL_REAL_F *y,
		     TSIL_REAL_F *z,
		     TSIL_REAL_F *u,
		     TSIL_REAL_F *v,
		     TSIL_REAL_F *qq,
		     TSIL_REAL_F *s)
{
  /* Below is cut and pasted directly from tsil_names.h */
  const char *uname[] = {"Uzxyv","Uuyxv","Uxzuv","Uyuzv"};
  const char *tname[] = {"Tvyz", "Tuxv", "Tyzv", "Txuv", "Tzyv", "Tvxu"};
  const char *sname[] = {"Svyz", "Suxv"};
  const char *bname[] = {"Bxz", "Byu"};
  const char *vname[] = {"Vzxyv","Vuyxv","Vxzuv","Vyuzv"};
  const char *tbarname[] = {"TBARvyz", "TBARuxv", "TBARyzv",
			    "TBARxuv", "TBARzyv", "TBARvxu"};

  TSIL_DATA res;
  int i, j;
  static int compiledSize = TSIL_REAL_SIZE_F;

  if (rsize_.rsize != compiledSize)
    TSIL_Warn ("tsilfevaluate", "Apparent type mismatch between Fortran and C!");

  TSIL_SetParameters (&res,
		      (TSIL_REAL) (*x),
		      (TSIL_REAL) (*y),
		      (TSIL_REAL) (*z), 
		      (TSIL_REAL) (*u),
		      (TSIL_REAL) (*v),
		      (TSIL_REAL) (*qq));

  TSIL_Evaluate (&res, (TSIL_REAL) (*s));

  /* Extract results */
  results_.M = (TSIL_COMPLEX_F) TSIL_GetFunction (&res, "M");

  /* Note non-standard indexing-by-pointer of UU, VV, TT, SS arrays
     for compatibility with Fortran indexing convention. */

  for (i=0; i<NUM_U_FUNCS; i++) {
    results_.U[i] = (TSIL_COMPLEX_F) TSIL_GetFunction (&res, uname[i]);
    results_.V[i] = (TSIL_COMPLEX_F) TSIL_GetFunction (&res, vname[i]);
    for (j=0; j<3; j++) {
      *(&(results_.UU[0][0]) + j * NUM_U_FUNCS + i) = 
        (TSIL_COMPLEX_F) TSIL_GetBoldFunction (&res, uname[i], j);
      *(&(results_.VV[0][0]) + j * NUM_U_FUNCS + i) = 
        (TSIL_COMPLEX_F) TSIL_GetBoldFunction (&res, vname[i], j);
    }
  }

  for (i=0; i<NUM_T_FUNCS; i++) {
    results_.T[i] = (TSIL_COMPLEX_F) TSIL_GetFunction (&res, tname[i]);
    results_.Tbar[i] = (TSIL_COMPLEX_F) TSIL_GetFunction (&res, tbarname[i]);
    for (j=0; j<3; j++)
      *(&(results_.TT[0][0]) + j * NUM_T_FUNCS + i) = 
        (TSIL_COMPLEX_F) TSIL_GetBoldFunction (&res, tname[i], j);
  }
  
  for (i=0; i<NUM_S_FUNCS; i++) {
    results_.S[i] = (TSIL_COMPLEX_F) TSIL_GetFunction (&res, sname[i]);
    for (j=0; j<3; j++)
      *(&(results_.SS[0][0])  + j * NUM_S_FUNCS + i) = 
        (TSIL_COMPLEX_F) TSIL_GetBoldFunction (&res, sname[i], j);
  }

  for (i=0; i<NUM_B_FUNCS; i++)
    results_.B[i] = (TSIL_COMPLEX_F) TSIL_GetFunction (&res, bname[i]);

  return;
}
