/* Analytic results for T-type functions. */

#include "internal.h"

/* ******************************************************************** */
/* A wrapper for the user API:                                          */

int TSIL_Tanalytic (TSIL_REAL X,
		    TSIL_REAL Y,   
		    TSIL_REAL Z,   
		    TSIL_COMPLEX S,
		    TSIL_REAL QQ,  
		    TSIL_COMPLEX *result)
{
  return Tanalytic (X, Y, Z, S, QQ, result);
}

/* ******************************************************************** */

int Tanalytic (TSIL_REAL X,
	       TSIL_REAL Y,
	       TSIL_REAL Z,
	       TSIL_COMPLEX S, 
	       TSIL_REAL QQ,
	       TSIL_COMPLEX *result)
{
  TSIL_REAL tmp;
  int success = 1;

  if (Y < Z) {tmp = Z; Z = Y; Y = tmp;}

  if (X < TSIL_TOL) {
    TSIL_Warn("Tanalytic", "T(x,y,z) is undefined for x = 0.");
    *result = TSIL_Infinity;
  }
  else if (Z < TSIL_TOL) *result = Tx0y (X,Y,S,QQ);
  else if (TSIL_CABS(S) < TSIL_TOL) *result = TAtZero(X, Y, Z, QQ);
  else if (TSIL_CABS(S-X) + TSIL_FABS(Y-Z) < TSIL_TOL) *result = TxyyAtx(X,Y,QQ);
  else if (TSIL_CABS(S-Z) + TSIL_FABS(Y-X) < TSIL_TOL) *result = TyyxAtx(Z,Y,QQ);
  else if (TSIL_CABS(S-Y) + TSIL_FABS(Z-X) < TSIL_TOL) *result = TyyxAtx(Y,X,QQ);
  else success = 0;

  return success;
}

/* ******************************************************************** */
/* hep-ph/0307101 eq. (6.15)                                            */

TSIL_COMPLEX Tx0y (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX sqDeltaSXY, tp, tm, log1mtp, log1mtm;
  TSIL_REAL    lnbarX, lnbarY;

  if (X < TSIL_TOL) { 
    TSIL_Warn("Tx0y", "T(x,0,y) is undefined for x = 0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(S) < TSIL_TOL) return -I2p(X, 0.0L, Y, QQ);
  if (Y < TSIL_TOL) return Tx00(X, S, QQ); 

  S = AddIeps(S);
  sqDeltaSXY = TSIL_CSQRT(X*X + Y*Y + S*S - 2.0L*X*Y - 2.0L*X*S - 2.0L*Y*S);
  lnbarX = TSIL_LOG(X/QQ); 
  lnbarY = TSIL_LOG(Y/QQ); 
  tp = (Y - X + S + sqDeltaSXY)/(2.0L * Y);
  tm = (Y - X + S - sqDeltaSXY)/(2.0L * Y);
  log1mtp = TSIL_CLOG(1.0L - tp);
  log1mtm = TSIL_CLOG(1.0L - tm);

  return (-Dilog (tp) - Dilog (tm) + (1.0L - Y/S)*log1mtp*log1mtm
	    + sqDeltaSXY * (log1mtp - log1mtm) / (2.0L * S)
	    + lnbarX*lnbarY - 0.5L*(lnbarY + 1.0L)*(lnbarY + 1.0L)
	    + (3.0L*S + Y - X) *(lnbarY-lnbarX)/(2.0L * S));
}


/* ******************************************************************** */
/* hep-ph/0307101 eq. (6.17)                                            */

TSIL_COMPLEX Tx00 (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_REAL lnbarX;  

  if (X < TSIL_TOL) {
    TSIL_Warn("Tx00", "T(x,0,0) is undefined for x = 0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(S) < TSIL_TOL) return -I2p(X, 0.0L, 0.0L, QQ);

  lnbarX = TSIL_LOG(X/QQ); 

  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL) 
    return 2.0L*Zeta2 -0.5L -lnbarX + 0.5L*lnbarX*lnbarX;

  S = AddIeps(S);

  return Zeta2 - 0.5L + (1.0L - X/S)*TSIL_CLOG(1.0L -S/X)
	    -lnbarX + 0.5L*lnbarX*lnbarX + Dilog (S/X);
}


/* ******************************************************************* */
/* SPM unpublished notes, probably equivalent to a derivative of
   F.A. Berends, A.I. Davydychev and N.I. Ussyukina,
   Phys. Lett. B426, 95 (1998) hep-ph/9712209,  eq. (25).             */
  
TSIL_COMPLEX TxyyAtx (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL QQ)
{
  TSIL_REAL lnbarX, lnbarY;
 
  if (X < TSIL_TOL) {
    TSIL_Warn("TxyyAtx", "T(x,y,y) is undefined for s = x = 0.");
    return TSIL_Infinity;
  }
  
  if (Y < TSIL_TOL) return Tx00(X,X,QQ);

  lnbarX = TSIL_LOG(X/QQ);
  lnbarY = TSIL_LOG(Y/QQ);   

  return (-0.5L + (Y/X-1.0L)*(Dilog(1.0 -X/Y) -Zeta2)
	    + lnbarX*(lnbarY - 1.0L) - 0.5L*lnbarY*lnbarY);
}


/* ******************************************************************* */
/* SPM unpublished notes, probably equivalent to a derivative of
   F.A. Berends, A.I. Davydychev and N.I. Ussyukina,
   Phys. Lett. B426, 95 (1998) hep-ph/9712209,  eq. (25).             */
    
TSIL_COMPLEX TyyxAtx (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL QQ)
{  
  TSIL_REAL lnbarX, lnbarY;

  if (X < TSIL_TOL) return -I2p(Y,Y,0,QQ);

  if (Y < TSIL_TOL) {
    TSIL_Warn("TyyxAtx", "T(y,y,x) with s = 0 is undefined for x = 0.");
    return TSIL_Infinity;
  }

  lnbarX = TSIL_LOG(X/QQ);
  lnbarY = TSIL_LOG(Y/QQ);
  return -0.5L + (1.0L - Y/X)*(Dilog(1.0 -X/Y) -Zeta2)
	 + lnbarX - 2.0L*lnbarY + 0.5L*lnbarY*lnbarY;
}

/* **************************************************************** */
                         
TSIL_COMPLEX TAtZero (TSIL_REAL x,
		      TSIL_REAL y,
		      TSIL_REAL z,
		      TSIL_REAL qq)
{
  /* If the first argument is 0, return TSIL_Infinity: */
  if (x < TSIL_TOL)
    return TSIL_Infinity;
  else
    return -I2p(x, y, z, qq);
} 
  
/* **************************************************************** */
  
TSIL_COMPLEX TprimeAtZero (TSIL_REAL x,
			   TSIL_REAL y,
			   TSIL_REAL z,
			   TSIL_REAL qq)
{
  /* If the first argument is 0, return TSIL_Infinity: */
  if (x < TSIL_TOL)
    return TSIL_Infinity;
  else
    return -(I2p2(x, y, z, qq) + x*I2p3(x, y, z, qq))/2.0L;
} 
