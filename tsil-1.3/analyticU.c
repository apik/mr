/* Analytic results for U-type functions */

#include "internal.h"

/* ******************************************************************* */

int TSIL_Uanalytic (TSIL_REAL X,
		    TSIL_REAL Y,
		    TSIL_REAL Z,
		    TSIL_REAL U, 
		    TSIL_COMPLEX S,
		    TSIL_REAL QQ,
		    TSIL_COMPLEX *result)
{
  int success = 1;

  if (TSIL_CABS(S) < TSIL_TOL) *result = TSIL_UAtZero (X, Y, Z, U, QQ);
  else if (X < TSIL_TOL) *result = TSIL_U0xyz(Y, Z, U, S, QQ);
  else if (U + TSIL_FABS(Y-Z) < TSIL_TOL) *result = TSIL_Uxy0y (X,Y,S,QQ);
  else if (Z + TSIL_FABS(Y-U) < TSIL_TOL) *result = TSIL_Uxy0y (X,Y,S,QQ);
  else if (Z + U < TSIL_TOL) *result = TSIL_Uxy00 (X,Y,S,QQ);
  else if (Y + U < TSIL_TOL) *result = TSIL_Ux00y (X,Z,S,QQ);
  else if (Y + Z < TSIL_TOL) *result = TSIL_Ux00y (X,U,S,QQ);
  else if (TSIL_CABS(S-X) + Y + TSIL_FABS(Z-U) < TSIL_TOL) *result = TSIL_Ux0yyAtx (X,Z,QQ);
  else if (TSIL_CABS(S-U) + Y + TSIL_FABS(X-Z) < TSIL_TOL) *result = TSIL_Uy0yxAtx (U,X,QQ);
  else if (TSIL_CABS(S-Z) + Y + TSIL_FABS(X-U) < TSIL_TOL) *result = TSIL_Uy0yxAtx (Z,X,QQ);
  else success = 0;

  return success;
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.25)                                           */

TSIL_COMPLEX TSIL_U0xyz (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL Z, 
			 TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_REAL temp, lnbarX, lnbarY, lnbarZ, DeltaXYZ;
  TSIL_COMPLEX log1mSoX, sqDeltaSYZ, sqDeltaXYZ;
  TSIL_COMPLEX Tp, Tm, logTp, logTm, log1mTp, log1mTm;
  TSIL_COMPLEX Rp, Rm, logRp, logRm, log1mRp, log1mRm, log1m1oRp, log1m1oRm;
  TSIL_COMPLEX part1, part2, part3;

  if (Y < Z) {temp=Z; Z=Y; Y=temp;}
 
  if (X < TSIL_TOL) return TSIL_U00xy (Y, Z, S, QQ);
  if (Z < TSIL_TOL) return TSIL_U0x0y (X, Y, S, QQ);
  if (TSIL_CABS(S) < TSIL_TOL) 
    return ((TSIL_I20xy(Y,Z,QQ) - TSIL_I2(X,Y,Z,QQ))/X);

  lnbarX = TSIL_LOG(X/QQ);
  lnbarY = TSIL_LOG(Y/QQ);
  lnbarZ = TSIL_LOG(Z/QQ);
  DeltaXYZ = X*X + Y*Y + Z*Z -2.0L*X*Y -2.0L*X*Z -2.0L*Y*Z;
  sqDeltaXYZ = TSIL_CSQRT(DeltaXYZ);
  Rp = (Y+Z-X + sqDeltaXYZ)/(2.0L*Y);
  Rm = (Y+Z-X - sqDeltaXYZ)/(2.0L*Y);
  logRp = TSIL_CLOG(Rp);
  logRm = TSIL_CLOG(Rm);
  log1mRp = TSIL_CLOG(1.0L - Rp);
  log1mRm = TSIL_CLOG(1.0L - Rm);

  /* DGR added in v1.2 to avoid branch problems if 1-1/Rp and/or
     1-1/Rm is real and negative */
  if (DeltaXYZ >= 0.0 && TSIL_CABS(Rp) < 1.0L)
    log1m1oRp = TSIL_CLOG (1.0L/Rp - 1.0L) + I*PI;
  else
    log1m1oRp = TSIL_CLOG (1.0L - 1.0L/Rp);

  if (DeltaXYZ >= 0.0 && TSIL_CABS(Rm) < 1.0L)
    log1m1oRm = TSIL_CLOG (1.0L/Rm - 1.0L) + I*PI;
  else
    log1m1oRm = TSIL_CLOG (1.0L - 1.0L/Rm);

  S = TSIL_AddIeps(S);
  sqDeltaSYZ = TSIL_CSQRT(S*S + Y*Y + Z*Z -2.0L*S*Y -2.0L*S*Z -2.0L*Y*Z);
  Tp = (S+Y-Z + sqDeltaSYZ)/(2.0L*Y); 
  Tm = (S+Y-Z - sqDeltaSYZ)/(2.0L*Y); 
  logTp = TSIL_CLOG(Tp);
  logTm = TSIL_CLOG(Tm);  
  log1mTp = TSIL_CLOG(1.0L - Tp);
  log1mTm = TSIL_CLOG(1.0L - Tm);

  part1 = (Z-Y)*(TSIL_Dilog(Tp) + TSIL_Dilog(Tm))/X -Z*log1mTp*log1mTm/S;

  if (TSIL_CABS(1-S/X) < 10.0L*TSIL_TOL)
    part2 = sqDeltaXYZ*(logRp - logRm)/(2.0L*X);
  else {
    log1mSoX = TSIL_CLOG(1.0L - S/X);

    part2 = (1.0L - X/S)*(
         ((Y-Z+sqDeltaXYZ)/(2.0L*X))*
	 (TSIL_Dilog(Rm) + TSIL_Dilog(1.0L/Rp) + TSIL_Dilog((Rp+Tp-1.0L)/Rp) 
	  + TSIL_Dilog((Rp+Tm-1.0L)/Rp)
	  -logRp*log1mSoX
	  - TSIL_CLOG(1.0L/Rm)*log1mRm - logRp*log1m1oRp
          - (log1mTp -logRp -TSIL_CLOG((1.0L-Tp)/Rp))*TSIL_CLOG((Rp+Tp-1.0L)/Rp)
	  - (log1mTm -logRp -TSIL_CLOG((1.0L-Tm)/Rp))*TSIL_CLOG((Rp+Tm-1.0L)/Rp))
	 + ((Y-Z-sqDeltaXYZ)/(2.0L*X))*(TSIL_Dilog(Rp) + TSIL_Dilog(1.0L/Rm)
          + TSIL_Dilog((Rm+Tp-1.0L)/Rm) + TSIL_Dilog((Rm+Tm-1.0L)/Rm) 
					- logRm*log1mSoX
         - TSIL_CLOG(1.0L/Rp)*log1mRp - logRm*log1m1oRm
          - (log1mTp -logRm -TSIL_CLOG((1.0L-Tp)/Rm))*TSIL_CLOG((Rm+Tp-1.0L)/Rm)
	  - (log1mTm -logRm -TSIL_CLOG((1.0L-Tm)/Rm))*TSIL_CLOG((Rm+Tm-1.0L)/Rm))
       + (0.5L*lnbarY + 0.5L*lnbarZ -2.0L)*log1mSoX 
       + 0.25L*log1mTp*log1mTp + 0.25L*log1mTm*log1mTm 
       + ( 0.5L*sqDeltaSYZ/(X-S) - 0.5L)*log1mTp  
       + (-0.5L*sqDeltaSYZ/(X-S) - 0.5L)*log1mTm  
	 + 2.0L*Zeta2*(Z-Y)/X );
  }

  part3 = -TSIL_I2(X,Y,Z,QQ)/X 
    + 0.5L*(lnbarY-lnbarZ)*(Z-Y+X +0.5L*X*(lnbarY-lnbarZ))/S +
    + 3.0L + lnbarY*(2.0L*Y/X - 1.0L) + 2.0L*Z*lnbarZ/X -2.5L*(Y+Z)/X
    + lnbarY*lnbarY*(0.5L*(Z-Y)/X - 0.25L) + lnbarY*lnbarZ*(0.5L - Z/X)   
    - 0.25L*lnbarZ*lnbarZ;

  return part1 + part2 + part3;
}

/* ******************************************************************* */
/* Special case of hep-ph/0307101 eq. (6.25)                           */

TSIL_COMPLEX TSIL_U00xy (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S,
			 TSIL_REAL QQ)
{
  TSIL_REAL lnbarX, lnbarY, temp;
  TSIL_COMPLEX lnbarmS, sqDeltaSXY, tp, tm, log1mtp, log1mtm;

  if (X < Y) {temp=Y; Y=X; X=temp;}

  if (X < TSIL_TOL) return TSIL_U0000(S, QQ);
  if (Y < TSIL_TOL) return TSIL_U000x(X, S, QQ);
  if (TSIL_FABS(1.0L-Y/X) < TSIL_TOL) return TSIL_U00xx(X, S, QQ);
  if (TSIL_CABS(S) < 10.0*TSIL_TOL){
    TSIL_Warn("TSIL_U00xy", "U(0,0,x,y) is undefined for s=0.");
    return TSIL_Infinity;
  }

  S = TSIL_AddIeps(S);
  lnbarX = TSIL_CLOG(X/QQ);
  lnbarY = TSIL_CLOG(Y/QQ);
  lnbarmS = TSIL_CLOG(-S/QQ);
  sqDeltaSXY = TSIL_CSQRT(S*S + X*X + Y*Y -2.0L*S*X -2.0L*S*Y -2.0L*X*Y);
  tp = (X-Y+S+sqDeltaSXY)/(2.0L*X);
  tm = (X-Y+S-sqDeltaSXY)/(2.0L*X);
  log1mtp = TSIL_CLOG(1.0L-tp);
  log1mtm = TSIL_CLOG(1.0L-tm);

  return (5.5L - lnbarmS + sqDeltaSXY*(log1mtm - log1mtp)/S 
          + (Y*(S-X+Y)*log1mtm*log1mtp)/(S*(X - Y)) 
          + (Y-X-S)*lnbarX/S + lnbarmS*lnbarX/(1.0L - Y/X)  
          + 0.5L*(X + Y)*lnbarX*lnbarX/(Y-X) 
          + (X-Y-S)*lnbarY/S + (lnbarmS-lnbarX)*lnbarY/(1.0L - X/Y) 
          + ((Y-X)/S + (X+Y)/(Y-X))*(TSIL_Dilog(tm)+ TSIL_Dilog(tp))   
          + (X-Y)*TSIL_Dilog(1.0L - Y/X)/S );
}


/* ******************************************************************* */
/* Special case of hep-ph/0307101 eq. (6.25)                           */

TSIL_COMPLEX TSIL_U0x0y (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S,
			 TSIL_REAL QQ)
{
  TSIL_COMPLEX lnbarXmS, result;
  TSIL_REAL lnbarY;

  if (X < TSIL_TOL)
    return TSIL_U000x(Y, S, QQ);

  if (Y < TSIL_TOL)
    return TSIL_U0x00(X, S, QQ);

  if (TSIL_CABS(S) < 10.0*TSIL_TOL)
    return ((TSIL_I200x(Y,QQ) - TSIL_I20xy(X,Y,QQ))/X);

  S = TSIL_AddIeps(S);    
  lnbarY = TSIL_LOG(Y/QQ);  
  result = 5.5L - 2.5L*(X + Y)/S - Zeta2*Y/X - TSIL_I20xy(X,Y,QQ)/S 
          - Y*TSIL_Dilog(S/Y)/X + Y*lnbarY*(1.0L/S - 0.5L*lnbarY/X) 
          + (Y/S - 1.0L)*TSIL_CLOG((Y - S)/QQ);

  if (TSIL_CABS(1 - S/X) > 10.0*TSIL_TOL) {
    lnbarXmS = TSIL_CLOG((X-S)/QQ);
    result += (1/X - 1/S)*((X-Y)*TSIL_Dilog((X-Y)/(X-S))
	      +lnbarXmS*(0.5L*(X-Y)*lnbarXmS + Y*lnbarY-2.0L*X));
  }
  return result;
}

/* ******************************************************************* */
/* Special case of hep-ph/0307101 eq. (6.25)                           */

TSIL_COMPLEX TSIL_U00xx (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX sqDeltaSXX, tp, log1mtp;
  TSIL_REAL lnbarX;

  if (X < TSIL_TOL)
    return TSIL_U0000(S, QQ);

  if (TSIL_CABS(S) < 10.0*TSIL_TOL){
    TSIL_Warn("TSIL_U00xx", "U(0,0,x,x) is undefined for s=0.");
    return TSIL_Infinity;
  }

  S = TSIL_AddIeps(S);  
  sqDeltaSXX = TSIL_CSQRT(S*S - 4.0L*X*S);
  lnbarX = TSIL_LOG(X/QQ);
  tp = 0.5L*(S+sqDeltaSXX)/X;
  log1mtp = TSIL_CLOG(1.0L - tp);

  return ((0.5L + X/S)*log1mtp - 3.0*sqDeltaSXX/S)*log1mtp 
          +5.5L -3.0L*lnbarX -0.5L*lnbarX*lnbarX + lnbarX*TSIL_CLOG(-S/QQ);
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.22)                                           */

TSIL_COMPLEX TSIL_Ux000 (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX log1mSoX;
  TSIL_REAL lnbarX;  

  if (X < TSIL_TOL)
    return TSIL_U0000 (S, QQ);

  lnbarX = TSIL_LOG(X/QQ);

  if (TSIL_CABS(S) < TSIL_TOL) 
    return (2.5L + Zeta2 - 2.0L*lnbarX + 0.5L*lnbarX*lnbarX);

  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL) 
    return (5.5L + 2.0L*Zeta2 - 3.0L*lnbarX + 0.5*lnbarX*lnbarX);

  S = TSIL_AddIeps(S);
  log1mSoX = TSIL_CLOG(1.0L - S/X);
  return (5.5L + Zeta2 + (1.0L-X/S)*log1mSoX*(log1mSoX - 3.0L + lnbarX)
	  -3.0L*lnbarX + 0.5L*lnbarX*lnbarX + TSIL_Dilog(S/X));
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.23)                                           */

TSIL_COMPLEX TSIL_U0x00 (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX log1mSoX, lnbarmS;
  TSIL_REAL lnbarX;  

  if (X < TSIL_TOL)
    return TSIL_U0000 (S, QQ);

  lnbarX = TSIL_LOG(X/QQ);

  if (TSIL_CABS(S) < TSIL_TOL) 
    return (5.5L + Zeta2  - 2.0L*lnbarX + 0.5L*lnbarX*lnbarX);
  
  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL) 
    return (5.5L + Zeta2 - 3.0L*lnbarX + 0.5L*lnbarX*lnbarX+I*PI);

  S = TSIL_AddIeps(S);
  log1mSoX = TSIL_CLOG(1.0L - S/X);
  lnbarmS = TSIL_CLOG(-S/QQ);

  return (5.5L + Zeta2 - lnbarmS -2.0L*lnbarX + 0.5L*lnbarX*lnbarX
	  + (1-X/S)*(log1mSoX*(lnbarmS-2.0L)+ TSIL_Dilog(S/X)));
}


/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.24)                                           */

TSIL_COMPLEX TSIL_U000x (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX log1mSoX;
  TSIL_REAL lnbarX;  

  if (X < TSIL_TOL)
    return TSIL_U0000 (S, QQ);

  if (TSIL_CABS(S) < TSIL_TOL){
    TSIL_Warn("TSIL_U000x", "U(0,0,0,x) is undefined for s=0.");
    return TSIL_Infinity;
  }

  lnbarX = TSIL_LOG(X/QQ);

  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL) 
    return (5.5L - 3.0L*Zeta2 + lnbarX*(0.5L*lnbarX - 3.0L - I*PI) +I*PI);

  S = TSIL_AddIeps(S);    
  log1mSoX = TSIL_CLOG(1.0L - S/X);

  return (5.5L - Zeta2 + (2.0L*X/S - 2.0L)*log1mSoX +
	  (lnbarX-1.0L)*TSIL_CLOG(-S/QQ) 
	  -2.0L*lnbarX - 0.5L*lnbarX*lnbarX -(1.0L +X/S)*TSIL_Dilog(S/X));
   
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.12)                                           */

TSIL_COMPLEX TSIL_U0000 (TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX lnbarmS;

  if (TSIL_CABS(S) < TSIL_TOL){
    TSIL_Warn("TSIL_U0000", "U(0,0,0,0) is undefined for s=0.");
    return TSIL_Infinity;
  }

  S = TSIL_AddIeps(S);
  lnbarmS = TSIL_CLOG(-S/QQ);
  return 0.5L*lnbarmS*lnbarmS - 3.0L*lnbarmS + 5.5L;
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.20)                                           */

TSIL_COMPLEX TSIL_Uxy0y (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S,
			 TSIL_REAL QQ)
{
  if (Y < TSIL_TOL) return TSIL_Ux000 (X, S, QQ);

  return (1.0L - TSIL_Tx0y(Y,X,S,QQ) 
	  + (2.0L - TSIL_LOG(Y/QQ))*TSIL_B(X,Y,S,QQ));
}

/* ********************************************************************** */
/* R. Scharf and J.B. Tausk, Nucl. Phys. B412, 523 (1994) eqs. (95),(83). */

TSIL_COMPLEX TSIL_Uxy00 (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S,
			 TSIL_REAL QQ)
{
  TSIL_COMPLEX sox, sqDeltasxy, r1, r2, onemr1, onemr2;
  TSIL_COMPLEX onem1or1, onem1or2, onemsox;
  TSIL_COMPLEX log1msox, logr1, logr2, lnbarx, lnbary, result;

  if (TSIL_FABS(X-Y) < TSIL_TOL) return TSIL_Uxx00 (X,S,QQ);
  if (X < TSIL_TOL) return TSIL_U0x00 (Y,S,QQ);  
  if (Y < TSIL_TOL) return TSIL_Ux000 (X,S,QQ);

  if (TSIL_CABS(S) < TSIL_TOL) 
    return (TSIL_I2(0,0,Y,QQ) - TSIL_I2(0,0,X,QQ))/(X-Y);

  S = TSIL_AddIeps(S);

  sox = S/X; /* Scharf and Tausk's x */
  sqDeltasxy = TSIL_CSQRT(S*S + X*X + Y*Y - 2.0L*(X*Y + X*S + Y*S));
  r1 = (X+Y-S + sqDeltasxy)/(2.0L*X);
  r2 = (X+Y-S - sqDeltasxy)/(2.0L*X);
           
  onemr1 = 1.0L - r1;
  onemr2 = 1.0L - r2;
  onem1or1 = 1.0L - 1/r1;
  onem1or2 = 1.0L - 1/r2;
  onemsox = 1.0L - sox;
  log1msox = TSIL_CLOG(onemsox);
  logr1 = TSIL_CLOG(r1);
  logr2 = TSIL_CLOG(r2);
  lnbarx = TSIL_CLOG(X/QQ);
  lnbary = TSIL_CLOG(Y/QQ);

  result = Zeta2 + 5.5L + (X/S -1.0L)*log1msox 
         + ((S + X - Y)/(4.0L*S))*lnbarx*lnbarx 
         + ((S - X + Y)/(4.0L*S))*lnbary*lnbary +
         + (X - Y - S)*lnbary/S 
         + (Y - X - 2.0L*S)*lnbarx/S 
         + (S + X - Y)/(2.0L*S)*TSIL_Dilog(sox) 
         + X*((r1 - r2)/(2.0L*S))*
           ((logr1 - logr2)*(2.0L - lnbary) 
           + TSIL_Dilog(onemr2) - TSIL_Dilog(onemr1) 
           + TSIL_Dilog(r1*onem1or2) - TSIL_Dilog(r2*onem1or1) 
           + TSIL_EtaBranch(onemsox, 1.0L/r2)*TSIL_CLOG(r1*onem1or2) 
           - TSIL_EtaBranch(onemsox, 1.0L/r1)*TSIL_CLOG(r2*onem1or1));

  return result;
}

/* ********************************************************************** */
/* R. Scharf and J.B. Tausk, Nucl. Phys. B412, 523 (1994) eqs. (94),(83). */

TSIL_COMPLEX TSIL_Uxx00 (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX sox, onemsox, sqrtonem4xos, r1, r2;
  TSIL_COMPLEX onemr1, onemr2, lnbarx;

  S = TSIL_AddIeps(S);

  if (TSIL_CABS(X/S) < TSIL_TOL) return TSIL_U0000(S, QQ);

  sox = S/X;
  lnbarx = TSIL_CLOG(X/QQ);

  if (TSIL_CABS(sox) < 10.0*TSIL_TOL) 
     return 0.5L + Zeta2 - lnbarx + 0.5L*lnbarx*lnbarx;

  onemsox = 1.0L - sox;

  sqrtonem4xos = TSIL_CSQRT(1.0L - 4.0L/sox);
  r1 = 1.0L + 0.5L*sox*(sqrtonem4xos -1.0L);
  r2 = 1.0L + 0.5L*sox*(-sqrtonem4xos -1.0L);
  onemr1 = 1.0L - r1;
  onemr2 = 1.0L - r2;

  return Zeta2 + 5.5L + 0.5L*lnbarx*lnbarx
           + lnbarx*(-3.0L - sqrtonem4xos*TSIL_CLOG(r1))
           + (1.0L/sox - 1.0L)*TSIL_CLOG(onemsox) + 0.5L*TSIL_Dilog(sox)
           + 0.5L*sqrtonem4xos*(
             4.0L*TSIL_CLOG(r1) - TSIL_Dilog(onemr1) + TSIL_Dilog(onemr2)
             - TSIL_Dilog(r2*onemr2) + TSIL_Dilog(r1*onemr1)
             - TSIL_EtaBranch(r2, onemsox)*TSIL_CLOG(r2*onemr2) +
             + TSIL_EtaBranch(r1, onemsox)*TSIL_CLOG(r1*onemr1));
}

/* ******************************************************************** */
/* hep-ph/0312092 eq. (A.17)                                            */

TSIL_COMPLEX TSIL_Ux00y (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX AX, AY, B0X, I200Y, S0XY, TX0Y, TY0X;

  if (Y < TSIL_TOL) return TSIL_Ux000(X, S, QQ);
  if (X < TSIL_TOL) return TSIL_U000x(Y, S, QQ);

  I200Y = TSIL_I2(0,0,Y,QQ);

  if (TSIL_CABS(S/X) < 10.0L*TSIL_TOL) return (I200Y - TSIL_I2(0,X,Y,QQ))/X; 

  AX = TSIL_A(X,QQ);
  AY = TSIL_A(Y,QQ);
  TX0Y = TSIL_Tx0y(X,Y,S,QQ);
  TY0X = TSIL_Tx0y(Y,X,S,QQ);

  if (TSIL_CABS(1.0L - S/X) < 10.0L*TSIL_TOL) return 
     (-2.0L*AY + 2.0L*AX*AY/X - TX0Y*Y - 2.0L*TY0X*Y)/Y;

  B0X = TSIL_B(0.0L,X,S,QQ);
  S0XY = TSIL_S0xy(X,Y,S,QQ);

  return (-AX - AY - I200Y + 2.0L*S0XY + 0.75L*X + 2.0L*TX0Y*X 
         + Y + TY0X*Y)/(S - X) -0.25L - TY0X - AY*B0X/Y;
}

/* ******************************************************************* */
/* SPM unpublished notes (is there a published source?)                */

TSIL_COMPLEX TSIL_Ux0yyAtx (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL QQ)
{
  TSIL_REAL lnbarX, lnbarY;
  
  if (X < 10.0*TSIL_TOL){
    TSIL_Warn("TSIL_U0x0yyAtx", "U(x,0,y,y) is undefined for s=x=0.");
    return TSIL_Infinity;
  }

  if (Y < TSIL_TOL) return TSIL_Ux000(X, X, QQ);

  lnbarX = TSIL_LOG(X/QQ);
  lnbarY = TSIL_LOG(Y/QQ);

  return (5.5L + (1.0L + Y/X)*Zeta2 - 3.0L*lnbarX + lnbarX*lnbarY
          - 0.5L*lnbarY*lnbarY - (1.0L +Y/X)*TSIL_Dilog(1.0L- X/Y)
          - 4.0L*TSIL_fPT(X,Y) );
}

/* ********************************************************************* */
/* SPM unpublished notes; this function might be useful for other things */

TSIL_COMPLEX TSIL_fPT(TSIL_REAL x, TSIL_REAL y)
{
  TSIL_REAL sqrtx, sqrty, argdilog;
 
  sqrtx = TSIL_SQRT(x);
  sqrty = TSIL_SQRT(y);
  argdilog = (sqrtx - sqrty)/(sqrtx + sqrty);

  return ((sqrty/sqrtx)*(1.5L*Zeta2 + TSIL_Dilog(argdilog) 
			 - TSIL_Dilog(-argdilog)));
}


/* ******************************************************************* */
/* SPM unpublished notes (is there a published source?)                */

TSIL_COMPLEX TSIL_Uy0yxAtx (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL QQ)
{
  TSIL_REAL lnbarX, lnbarY;
  TSIL_COMPLEX  log1mXoY;

  if (X < TSIL_TOL) return TSIL_Ux00y(Y, Y, 0, QQ);
  if (Y < TSIL_TOL) return TSIL_U000x(X, X, QQ);

  lnbarX = TSIL_LOG(X/QQ);
  lnbarY = TSIL_LOG(Y/QQ);
  log1mXoY = TSIL_CLOG(1.0L - TSIL_AddIeps(X/Y));

  return 5.5L - 2.0L*lnbarX - lnbarY + lnbarY*lnbarY/2.0L
          + TSIL_Dilog(1-X/Y)*(1.0L + 2.0L*Y/X) -2.0L*(1 + Y/X)*Zeta2
          + log1mXoY*(lnbarX - 1.0L + Y*(1.0L - lnbarY)/X);
}

/* ******************************************************************* */
/* Corrects the U values in cases where there were uncanceled double   */
/* poles in their derivatives.                                         */

void TSIL_CorrectUs (TSIL_DATA *foo)
{
  TSIL_REAL x, y, z, s, qq;
  int i;
  
  /* For convenience */
  s  = foo->s;
  qq = foo->qq;

  /* Loop over the Us, checking whether arg[1] is zero: */
  for (i=0; i<4; i++) {
    if (foo->U[i].arg[1] < TSIL_TOL) {
      x = foo->U[i].arg[0];
      y = foo->U[i].arg[2];
      z = foo->U[i].arg[3];

      if (TSIL_FABS(s-x) > 0.0L) {
	/* Generic case: uses eqns A17 and A18 from hep-ph/0312092 */

	if (TSIL_FABS(y - z) > 0.0L) { /* Use A17 */
	  foo->U[i].value = 
	    (1.0L/(z-y) + 1.0L/(s-x)) * y * *(foo->U[i].tval[2]) +
	    (1.0L/(y-z) + 1.0L/(s-x)) * z * *(foo->U[i].tval[1]) +
	    (2.0L * x * *(foo->U[i].tval[0]) + 2.0L * *(foo->U[i].sval) - 
	     TSIL_I2(0, y, z, qq) - TSIL_A(x, qq) - TSIL_A(y, qq) - TSIL_A(z, qq) 
             + x + y + z - s/4.0L)/(s - x) +
	     TSIL_B0x (x, s, qq)*(TSIL_A(y, qq) - TSIL_A(z, qq))/(z - y);
	}
	else { /* Use A18 */
	  foo->U[i].value =
	    ((s - x - 4.0L*y) * *(foo->U[i].tval[2]) -
	     4.0L * x * *(foo->U[i].tval[0]) -
	     4.0L * *(foo->U[i].sval) +
	     2.0L*TSIL_I2(0, y, y, qq) + 2.0L*TSIL_A(x, qq) + 4.0L*TSIL_A(y, qq) 
             - 3.0L*x - 4.0L*y + 1.5L*s)/(x - s) 
	     - (1.0L + TSIL_A(y, qq)/y)*TSIL_B0x (x, s, qq);
	}
      }
      else {
      /* The s = x case */
        if (TSIL_FABS(y - z) > 0.0L)
  	  foo->U[i].value = TSIL_Ux0yzAtx (x, y, z, foo->U[i].tval, qq);
        else
	  foo->U[i].value = TSIL_Ux0yyAtx (x, y, qq );    
      }
    }
  }

  return;
}

/* ******************************************************************* */
/* SPM unpublished notes                                               */

TSIL_COMPLEX TSIL_Ux0yzAtx (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_REAL z,
			    TSIL_COMPLEX **tptr,
			    TSIL_REAL qq)
{
  return (-2.0L*TSIL_A(y, qq) + 2.0L*TSIL_A(x, qq)*TSIL_A(y, qq)/x + 
	  2.0L*TSIL_A(z, qq) - 2.0L*TSIL_A(x, qq)*TSIL_A(z, qq)/x +
	  (z - y) * *(tptr[0]) - 2.0L * y * *(tptr[2]) 
          + 2.0L * z * *(tptr[1]))/(y - z);
}

/* **************************************************************** */
                         
TSIL_COMPLEX TSIL_UAtZero (TSIL_REAL x,   
			   TSIL_REAL y,
			   TSIL_REAL z,
			   TSIL_REAL u,
			   TSIL_REAL qq)
{
  if (TSIL_FABS(x - y) > TSIL_TOL)
    return (TSIL_I2(x, z, u, qq) - TSIL_I2(y, z, u, qq))/(y - x);
  else
    return -TSIL_I2p(x, z, u, qq);
}
    
/* **************************************************************** */
      
TSIL_COMPLEX TSIL_UprimeAtZero (TSIL_REAL x,
				TSIL_REAL y,
				TSIL_REAL z,
				TSIL_REAL u,
				TSIL_REAL qq)
{
  TSIL_REAL tmp;

  if (x < TSIL_TOL) {
    if (z < u) {tmp = u; u = z; z = tmp;}
    return z*TSIL_I2p2(z,0,u,qq);
  }  
    
  if (TSIL_FABS(x - y) > TSIL_TOL)
    return x*(TSIL_I2(x, z, u, qq) - TSIL_I2(y, z, u, qq))/TSIL_POW(y - x, 3)
      + x*TSIL_I2p(x, z, u, qq)/TSIL_POW(y - x, 2)
      + x*TSIL_I2p2(x, z, u, qq)/(2.0L*(y - x));
  else
    return -x*TSIL_I2p3(x, z, u, qq)/6.0L;
}
