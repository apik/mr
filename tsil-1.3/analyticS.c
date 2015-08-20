/* Analytic results for S-type functions. */

#include "internal.h"

/* ****************************************************************** */

int TSIL_Sanalytic (TSIL_REAL X,
		    TSIL_REAL Y,
		    TSIL_REAL Z,
		    TSIL_COMPLEX S,
		    TSIL_REAL QQ,
		    TSIL_COMPLEX *result)
{
  TSIL_REAL tmp;
  int success = 1;

  if (Z > Y) {tmp = Y; Y = Z; Z = tmp;}
  if (Z > X) {tmp = X; X = Z; Z = tmp;}

  if (TSIL_CABS(S) < TSIL_TOL)
    *result = TSIL_SAtZero (X,Y,Z,QQ);
  else if (Z < TSIL_TOL)
    *result = TSIL_S0xy (X, Y, S, QQ);
  else if (TSIL_CABS(S-X) +TSIL_FABS(Y-Z) < TSIL_TOL)
    *result = TSIL_SxyyAtx (X,Y,QQ);
  else if (TSIL_CABS(S-Y) +TSIL_FABS(X-Z) < TSIL_TOL)
    *result = TSIL_SxyyAtx (Y,X,QQ);
  else if (TSIL_CABS(S-Z) +TSIL_FABS(X-Y) < TSIL_TOL)
    *result = TSIL_SxyyAtx (Z,X,QQ);
  else success = 0;

  return success;
}

/* ****************************************************************** */
/* hep-ph/0307101 eq. (6.14)                                          */

TSIL_COMPLEX TSIL_S0xy (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S,
			TSIL_REAL QQ)
{
  TSIL_COMPLEX sqDeltaSXY, tp, tm, log1mtp, log1mtm;
  TSIL_REAL    lnbarX, lnbarY, temp;

  if (TSIL_CABS(S) < TSIL_TOL) return TSIL_I20xy(X, Y, QQ);

  if (X < Y) {temp = Y; Y = X; X = temp;}

  if (X < TSIL_TOL) return TSIL_S000 (S, QQ);
  if (Y < TSIL_TOL) return TSIL_S00x (X, S, QQ);

  S = TSIL_AddIeps(S);
  sqDeltaSXY = TSIL_CSQRT(X*X + Y*Y + S*S - 2.0L*X*Y - 2.0L*X*S - 2.0L*Y*S);
  lnbarX = TSIL_LOG(X/QQ);
  lnbarY = TSIL_LOG(Y/QQ);
  tp = (X - Y + S + sqDeltaSXY)/(2.0L * X);
  tm = (X - Y + S - sqDeltaSXY)/(2.0L * X);
  log1mtp = TSIL_CLOG(1.0L - tp);
  log1mtm = TSIL_CLOG(1.0L - tm);
  
  return ((Y-X)*(TSIL_Dilog(tp) + TSIL_Dilog(tm))
	  - Y*(1.0L-X/S)*log1mtp*log1mtm 
	  - 0.25L*(S+X+Y)*sqDeltaSXY*(log1mtp - log1mtm)/S
	  + 0.5L*(Y-X)*lnbarX*lnbarX - Y*lnbarX*lnbarY
	  + 0.25L*(Y*Y - X*X)*(lnbarX - lnbarY)/S
	  + (2.0L*X -0.25L*S)*lnbarX + (2.0L*Y -0.25L*S)*lnbarY
	  + 1.625L*S - 2.0L*X -2.0L*Y);
}  

/* ****************************************************************** */
/* hep-ph/0307101 eq. (6.16)                                          */

TSIL_COMPLEX TSIL_S00x (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_REAL lnbarX;  

  if (TSIL_CABS(S) < TSIL_TOL) return TSIL_I200x(X, QQ);
  if (X < TSIL_TOL) return TSIL_S000 (S, QQ);

  lnbarX = TSIL_LOG(X/QQ);

  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL)
    return (X*(-0.375L -2.0L*Zeta2 + (1.5L -0.5L*lnbarX)*lnbarX));

  S = TSIL_AddIeps(S);
  return (1.625L*S - (2.0L + Zeta2 + TSIL_Dilog(S/X))*X  
	  + 0.5L*(X*X/S - S)*TSIL_CLOG(1-S/X) 
	  + (2.0L*X - 0.5L*S)*lnbarX - 0.5L*X*lnbarX*lnbarX);
}

/* ****************************************************************** */
/* hep-ph/0307101 eq. (6.10)                                          */

TSIL_COMPLEX TSIL_S000 (TSIL_COMPLEX S, TSIL_REAL QQ)
{
  if (TSIL_CABS(S) < TSIL_TOL) return 0.0L + I*0.0L;

  S = TSIL_AddIeps(S);
  return (S*(1.625L - 0.5L*TSIL_CLOG(-S/QQ)));
  
}

/* ****************************************************************** */
/* SPM unpublished notes, probably equivalent to
   F.A. Berends, A.I. Davydychev and N.I. Ussyukina,
   Phys. Lett. B426, 95 (1998) hep-ph/9712209,  eq. (25).             */
  
TSIL_COMPLEX TSIL_SxyyAtx (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL QQ)
{ 
  TSIL_REAL lnbarX, lnbarY;

  if (X < TSIL_TOL) return TSIL_I2(0,Y,Y,QQ);
  if (Y < TSIL_TOL) return TSIL_S00x(X,X,QQ);

  lnbarX = TSIL_LOG(X/QQ);
  lnbarY = TSIL_LOG(Y/QQ);
  return -0.375L*X -4.0L*Y + (0.5L*X - Y)*lnbarY*lnbarY
         + (X + Y*Y/X - 2.0L*Y)*(TSIL_Dilog(1.0 -X/Y) -Zeta2)
         + (1.5L*X - Y)*lnbarX + 5.0L*Y*lnbarY - X*lnbarX*lnbarY;

}

/* ***************************************************************** */
                 
TSIL_COMPLEX TSIL_SAtZero (TSIL_REAL x,
			   TSIL_REAL y,
			   TSIL_REAL z,
			   TSIL_REAL qq)
{
  return TSIL_I2(x, y, z, qq);
}
  
/* ***************************************************************** */
  
TSIL_COMPLEX TSIL_SprimeAtZero (TSIL_REAL x,
				TSIL_REAL y,
				TSIL_REAL z,
				TSIL_REAL qq)
{
  TSIL_REAL tmp; 

  if (x < y) {tmp = y; y = x; x = tmp;}
  if (x < z) {tmp = z; z = x; x = tmp;}

  return (x/2.0L)*TSIL_I2p2(x, y, z, qq) - 0.125L;
}
