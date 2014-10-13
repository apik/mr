/* Definitions needed by the Fortran interface (fevaluate.c) */

#ifndef TSIL_FORTRAN_H
#define TSIL_FORTRAN_H

/*
   These values must match the types used in the calling Fortran
   programs.  TSIL_REAL_SIZE_F should be set to the size in bytes of
   the TSIL_REAL_F type.  The calling Fortran program should then use
   REAL type

   REAL*(TSIL_FORT_SIZE) 
*/

/* YOU MAY WISH TO UNCOMMENT THE FOLLOWING IF YOUR FORTRAN COMPILER
   SUPPORTS REAL*16/COMPLEX*32 */
/*  #if defined(TSIL_SIZE_LONG) */
/*  #define TSIL_REAL_F       long double */
/*  #define TSIL_COMPLEX_F    long double complex */
/*  #define TSIL_REAL_SIZE_F  16 */

/* Sets REAL*8 as the type used by the Fortran program: */
#define TSIL_REAL_F       double
#define TSIL_COMPLEX_F    double complex
#define TSIL_REAL_SIZE_F  8

#endif /* tsil_fortran.h */
