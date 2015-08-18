/* Analytic results for Tbar type functions */

#include "internal.h"

/* ************************************************************** */
/* This function is defined in hep-ph/0307101 eq. (2.21)          */

int TSIL_Tbaranalytic (TSIL_REAL X,
		       TSIL_REAL Y,
		       TSIL_REAL Z, 
		       TSIL_COMPLEX S,
		       TSIL_REAL QQ,
		       TSIL_COMPLEX *result)
{
  TSIL_REAL tmp;
  int success = 1;
  
  if (Y < Z) {tmp = Z; Z = Y; Y = tmp;}

  if (X < TSIL_TOL) 
    *result = TSIL_Tbar0xy(Y,Z,S,QQ);
  else if (Z < TSIL_TOL) 
    *result = TSIL_Tx0y(X,Y,S,QQ) + TSIL_LOG(X/QQ)*TSIL_B(0,Y,S,QQ);
  else 
    success = 0;

  return success;
}

/* ************************************************************** */
/* hep-ph/0307101 eq. (6.18)                                      */

TSIL_COMPLEX TSIL_Tbar0xy (TSIL_REAL X, TSIL_REAL Y, 
			   TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX sqDeltaSXY, tp, tm, log1mtp, log1mtq;
  TSIL_REAL    AX, AY, lnbarX, lnbarY, Sthresh, lnbarS, logXoverY, temp;
  TSIL_REAL    onemYoX, onemYoX2, onemYoX3, onemYoX4;

  if (X < Y) {temp = Y; Y = X; X = temp;}

  if (X < TSIL_TOL) return TSIL_Tbar000 (S, QQ);
  if (Y < TSIL_TOL) return TSIL_Tbar00x (X, S, QQ);

  if (TSIL_CABS(S) < TSIL_TOL) {
    onemYoX = 1.0L - Y/X;
    lnbarX = TSIL_LOG(X/QQ);

    if (onemYoX > 0.01) {
      AX = X*lnbarX-X;
      AY = Y*TSIL_LOG(Y/QQ)-Y;
      return ((X+Y)*TSIL_I20xy(X,Y,QQ) + 2.0L*(AX-Y)*(AY-X) + 
	      X*X + Y*Y)/((X-Y)*(X-Y));
    } 
    else {
      onemYoX2 = onemYoX*onemYoX; 
      onemYoX3 = onemYoX*onemYoX2; 
      onemYoX4 = onemYoX*onemYoX3; 

      return (-1.5L -lnbarX - 0.5L*lnbarX*lnbarX + 
	      (0.5L + 0.5L*lnbarX)*onemYoX + 
	      (0.02777777777777777777777777777777L + 
	       0.16666666666666666666666666666667L*lnbarX)*onemYoX2 + 
	      (-0.0138888888888888888888888888889L + 
	       0.0833333333333333333333333333333L*lnbarX)*onemYoX3 + 
	      (-0.0175L + 0.05L*lnbarX)*onemYoX4 + 
	      (-0.0155555555555555555555555555556L + 
	       0.0333333333333333333333333333333L*lnbarX)*onemYoX4*onemYoX + 
	      (-0.0130385487528344671201814058957L + 
	       0.0238095238095238095238095238092L*lnbarX)*onemYoX4*onemYoX2 + 
	      (-0.0108418367346938775510204081633L + 
	       0.0178571428571428571428571428571L*lnbarX)*onemYoX4*onemYoX3 + 
	      (-0.0090663580246913580246913580247L + 
	       0.0138888888888888888888888888889L*lnbarX)*onemYoX4*onemYoX4);
    }
  }

  lnbarX = TSIL_LOG(X/QQ);
  lnbarY = TSIL_LOG(Y/QQ);   

  if (TSIL_CABS(1.0L - S/(X+Y+2.0L*TSIL_SQRT(X*Y))) < TSIL_TOL) {
    Sthresh = X+Y+2.0L*TSIL_SQRT(X*Y);
    lnbarS = TSIL_LOG(Sthresh/QQ);

    return (2.0L*(Y-X)*TSIL_Dilog(TSIL_CSQRT(Y/Sthresh))/Sthresh +
	    (-Sthresh + 0.5L*lnbarX*lnbarX*(Y-Sthresh) 
	     + 0.5L*lnbarY*lnbarY*(X-Sthresh) 
	     + 0.5L*lnbarS*lnbarS*(Y-X) + lnbarX*(Sthresh+X-Y) 
	     + lnbarY*(Y+Sthresh-X) + lnbarS*lnbarY*(X-Y) - lnbarX*lnbarY*X  
	     + (8.L*X + 4.L*Y - 12.L*Sthresh)*Zeta2)/(2.0L*Sthresh));        
  }

  S = TSIL_AddIeps(S);
  sqDeltaSXY = TSIL_CSQRT(TSIL_Delta (S, X, Y));

  tm = (X - Y + S - sqDeltaSXY) / (2.0L * X);
  tp = (X - Y + S + sqDeltaSXY) / (2.0L * X);
  log1mtp = TSIL_CLOG(1.0L - tp);
  log1mtq = TSIL_CLOG((1.0L - tm)*X/sqDeltaSXY);
  logXoverY = lnbarX - lnbarY;

  return ( ( (Y-X-sqDeltaSXY)*TSIL_Dilog(tp)
	     + (Y-X+sqDeltaSXY)*TSIL_Dilog(tm)
	     + 2.0L*sqDeltaSXY*TSIL_Dilog(X*(tp-1.0L)/sqDeltaSXY)
	     + (S-X +2.0L*sqDeltaSXY)*log1mtp*log1mtp
	     + sqDeltaSXY*log1mtq*log1mtq
	     + 4.0L*sqDeltaSXY*log1mtq*log1mtp
	     + (sqDeltaSXY*(1.0L -lnbarY) + (S-X)*logXoverY)*log1mtp
	     - 2.0L*sqDeltaSXY*logXoverY*TSIL_CLOG(sqDeltaSXY/X)
	     + 0.5L*(lnbarX-1.0L)*(S*(1.0L-lnbarY) + (Y-X)*logXoverY)
	     + sqDeltaSXY*(2.0L*Zeta2 -1.5L*lnbarX*lnbarX +0.5L*lnbarX
		  + 2.5L*lnbarX*lnbarY -lnbarY*lnbarY - 0.5L*lnbarY))/S);
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.10)                                           */

TSIL_COMPLEX TSIL_Tbar000 (TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX lnbarmS;

  if (TSIL_CABS(S) < TSIL_TOL) {
     TSIL_Warn("TSIL_Tbar000", "Tbar(0,0,0) is undefined at s=0."); 
     return TSIL_Infinity;
  }

  S = TSIL_AddIeps(S);
  lnbarmS = TSIL_CLOG(-S/QQ);
  return -0.5L*lnbarmS*lnbarmS + lnbarmS - 0.5L; 
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.19)                                           */

TSIL_COMPLEX TSIL_Tbar00x (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX log1mSoX;
  TSIL_REAL    lnbarX;  

  if (X < TSIL_TOL) return TSIL_Tbar000 (S, QQ);

  lnbarX = TSIL_LOG(X/QQ); 

  if (TSIL_CABS(S) < TSIL_TOL) 
    return (0.5L - Zeta2 -0.5L*lnbarX*lnbarX);

  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL) 
    return (-0.5L -2.0L*Zeta2 + lnbarX - 0.5L*lnbarX*lnbarX);

  S = TSIL_AddIeps(S);
  log1mSoX = TSIL_CLOG(1.0L - S/X);
  return (-TSIL_Dilog(S/X) -0.5L - Zeta2 + lnbarX - 0.5L*lnbarX*lnbarX
	  + (1.0L - X/S)*log1mSoX*(1.0L -lnbarX - log1mSoX));
}
