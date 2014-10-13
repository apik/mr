/* Analytic results for M */

#include "internal.h"

/* ******************************************************************* */
/* A wrapper for the user API:                                         */

int TSIL_Manalytic (TSIL_REAL X,
		    TSIL_REAL Y,
		    TSIL_REAL Z,
		    TSIL_REAL U,
		    TSIL_REAL V,
		    TSIL_COMPLEX S,
		    TSIL_COMPLEX *result)
{
  return Manalytic (X, Y, Z, U, V, S, result);
}

/* ******************************************************************* */

int Manalytic (TSIL_REAL X,
	       TSIL_REAL Y,
	       TSIL_REAL Z,
	       TSIL_REAL U,
	       TSIL_REAL V,
               TSIL_COMPLEX S,
	       TSIL_COMPLEX *result)
{
  int success = 1;

  if (V + TSIL_FABS(X-Y) + TSIL_FABS(Z-U) < TSIL_TOL) *result = Mxxyy0(X, Z, S);
  else if (X + Y + Z + U < TSIL_TOL) *result = M0000x(V, S);
  else if (X + Y + V < TSIL_TOL) *result = M00xy0(Z, U, S);
  else if (Z + U + V < TSIL_TOL) *result = M00xy0(X, Y, S);
  else if (X + Z + TSIL_FABS(Y-V) + TSIL_FABS(U-V) < TSIL_TOL) *result = M0x0xx(Y, S);
  else if (Y + U + TSIL_FABS(X-V) + TSIL_FABS(Z-V) < TSIL_TOL) *result = M0x0xx(X, S);
  else if (X + Z + V < TSIL_TOL) *result = M0x0y0(Y, U, S);
  else if (Y + U + V < TSIL_TOL) *result = M0x0y0(X, Z, S);
  else if (X + Y + Z + TSIL_FABS(U-V) < TSIL_TOL) *result = M000xx(U, S);
  else if (X + Y + U + TSIL_FABS(Z-V) < TSIL_TOL) *result = M000xx(Z, S);
  else if (X + Z + U + TSIL_FABS(Y-V) < TSIL_TOL) *result = M000xx(Y, S);
  else if (Y + Z + U + TSIL_FABS(X-V) < TSIL_TOL) *result = M000xx(X, S);
  else if (X + U + V + TSIL_FABS(Y-Z) < TSIL_TOL) *result = M0xx00(Y, S);
  else if (Y + Z + V + TSIL_FABS(X-U) < TSIL_TOL) *result = M0xx00(X, S);
  else if (X + V + TSIL_FABS(Y-U) + TSIL_FABS(Z-U) < TSIL_TOL) *result = M0xxx0(Y,S);
  else if (Y + V + TSIL_FABS(X-U) + TSIL_FABS(Z-U) < TSIL_TOL) *result = M0xxx0(X,S);
  else if (Z + V + TSIL_FABS(X-U) + TSIL_FABS(Y-U) < TSIL_TOL) *result = M0xxx0(X,S);
  else if (U + V + TSIL_FABS(X-Z) + TSIL_FABS(Y-Z) < TSIL_TOL) *result = M0xxx0(X,S);
  else if (TSIL_CABS(S-V) + X + U + TSIL_FABS(Z-Y) < TSIL_TOL) *result = M0yy0xAtx(V,Y);
  else if (TSIL_CABS(S-V) + Y + Z + TSIL_FABS(X-U) < TSIL_TOL) *result = M0yy0xAtx(V,X);
  else if (TSIL_CABS(S-Z) + X + U + TSIL_FABS(Y-V) < TSIL_TOL) *result = M0xy0yAtx(Z,Y);
  else if (TSIL_CABS(S-X) + Z + Y + TSIL_FABS(U-V) < TSIL_TOL) *result = M0xy0yAtx(X,U);
  else if (TSIL_CABS(S-Y) + U + X + TSIL_FABS(Z-V) < TSIL_TOL) *result = M0xy0yAtx(Y,Z);
  else if (TSIL_CABS(S-U) + Y + Z + TSIL_FABS(X-V) < TSIL_TOL) *result = M0xy0yAtx(U,X);
  else if (TSIL_CABS(S) < TSIL_TOL) *result = sMprimeAtZero(X, Y, Z, U, V);
  else success = 0;

  return success;
}

/* ******************************************************************* */
/* Defined in D.J. Broadhurst, Z. Phys. C47, 115 (1990), eq. (13).     */

TSIL_COMPLEX F3plusBroadhurst (TSIL_COMPLEX z)
{
  TSIL_COMPLEX logz = TSIL_CLOG (z);
  return (6.0L*Trilog(z) - 4.0L*logz*Dilog(z) - TSIL_CLOG (1.0L - z)*logz*logz);
}

/* ******************************************************************* */
/* Defined in D.J. Broadhurst, Z. Phys. C47, 115 (1990), eq. (13).     */

TSIL_COMPLEX F3minusBroadhurst (TSIL_COMPLEX z)
{
  TSIL_COMPLEX logz = TSIL_CLOG (z);
  return (6.0L*Trilog(-z) - 4.0L*logz*Dilog(-z) - TSIL_CLOG (1.0L + z)*logz*logz);
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.13)                                           */

TSIL_COMPLEX M00000 (TSIL_COMPLEX S)
{
  if (TSIL_CABS(S) < TSIL_TOL) {
    TSIL_Warn("M00000", "M(0,0,0,0,0) is undefined for s=0."); 
    return TSIL_Infinity;
  }

  return -6.0L*Zeta3/S;
}

/* ******************************************************************* */
/* D.J. Broadhurst, Z. Phys. C47, 115 (1990), eq. (20).                */

TSIL_COMPLEX M00xx0 (TSIL_REAL X, TSIL_COMPLEX S)
{
  if (X < TSIL_TOL) return M00000(S);
 
  if (TSIL_CABS(S) < TSIL_TOL) return ((I2(0.,X,X,X) - 2.0L*I2(0.,0.,X,X))/(X*X));

  if (TSIL_CABS(1.0L - S/X) < 10.0*TSIL_TOL) {
    TSIL_Warn("M00xx0", "M(0,0,x,x,0) is undefined for s=x."); 
    return TSIL_Infinity;
  }

  S = AddIeps(S);
  return (F3plusBroadhurst(X/(X-S)) - 6.0L*Zeta3)/S;
}

/* ******************************************************************* */
/* D.J. Broadhurst, Z. Phys. C47, 115 (1990), eq. (14).                */

TSIL_COMPLEX M000x0 (TSIL_REAL X, TSIL_COMPLEX S)
{
  if (X < TSIL_TOL) return M00000(S);
 
  if (TSIL_CABS(S) < TSIL_TOL) {
    TSIL_Warn("M000x0", "M(0,0,0,x,0) is undefined for s=0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(1.0L - S/X) < 10.0*TSIL_TOL) 
    return ((-3.0L*Zeta3 + I*PI*PI*PI/3.0L)/X);

  S = AddIeps(S);
  return (F3plusBroadhurst(X/(X-S)) - F3plusBroadhurst(S/(S-X)) -
	   6.0L*Zeta3)/(2.0L*S);
}

/* ******************************************************************* */
/* R. Scharf and J.B. Tausk, Nucl. Phys. B412, 523 (1994). eq. (104)   */

TSIL_COMPLEX M0000x (TSIL_REAL X, TSIL_COMPLEX S)
{
  TSIL_COMPLEX SoX, SoXp1, log1mSoX, logmSoX;

  if (X < TSIL_TOL) return M00000(S);
 
  if (TSIL_CABS(S) < TSIL_TOL) {
    TSIL_Warn("M0000x", "M(0,0,0,0,x) is undefined for s=0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(1.0L - S/X) < 10.0*TSIL_TOL) return ((-3.0L*Zeta3 + I*PI*Zeta2)/S);

  S = AddIeps(S);
  SoX = S/X;
  SoXp1 = SoX + 1.0L;
  log1mSoX = TSIL_CLOG(1.0L - SoX);
  logmSoX = TSIL_CLOG(-SoX);

  return ((2.0L*Zeta2*(logmSoX - log1mSoX) 
	   + (2.0L/3.0L)*log1mSoX*log1mSoX*log1mSoX
	   - TSIL_CLOG(SoXp1)*logmSoX*logmSoX 
	   + (4.0L*log1mSoX - 2.0L*logmSoX)*Dilog(SoXp1)
	   - 4.0L*Trilog(SoXp1) - 2.0L*Trilog(-SoX)
	   - 2.0L*Trilog(SoX) - 4.0L*Trilog(X/(X-S))
	   - 4.0L*Trilog((S+X)/(S-X)) + 4.0L*Trilog((S+X)/(X-S)) + Zeta3)/S);
}

/* ******************************************************************* */
/* J. Fleischer, A.V. Kotikov and O.L. Veretin,                        */
/* Nucl. Phys. B547, 343 (1999) [hep-ph/9808242] I_125, eq. (79)       */

TSIL_COMPLEX M0x0xx (TSIL_REAL X, TSIL_COMPLEX S)
{
  TSIL_COMPLEX sqrt1m4XoS, sm1osp1, logsm1osp1;

  if (X < TSIL_TOL) return M00000(S);
 
  if (TSIL_CABS(S) < TSIL_TOL) {
    TSIL_Warn("M0x0xx", "M(0,x,0,x,x) is undefined for s=0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(1.0L - 0.25L*S/X) < 10.0*TSIL_TOL)
    return ((-2.625L*Zeta3 + 3.0L*Zeta2*TSIL_LOG(2.0L) + I*PI*PI*PI/8.0L)/X);

  S = AddIeps(S);
  sqrt1m4XoS = TSIL_CSQRT(1.0L - 4.0L*X/S);
  sm1osp1 = (sqrt1m4XoS - 1.0L)/(sqrt1m4XoS +1.0L);
  logsm1osp1 = TSIL_CLOG(sm1osp1);

  return (6.0L*Trilog(sm1osp1) - 6.0L*logsm1osp1*Dilog(sm1osp1)
	  -2.0*logsm1osp1*logsm1osp1*TSIL_CLOG(1.0L-sm1osp1) - 6.0L*Zeta3)/S;
}

/* ******************************************************************* */
/* Don't know a published source for this one. Maybe an SPM original?  */

TSIL_COMPLEX M00xy0 (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S)
{
  TSIL_COMPLEX YoX, log1mYoX, XoY, logXoY, onemSoX, onemSoY, SmYoSmX;
  TSIL_COMPLEX log1mSoX, log1mSoY;
  TSIL_REAL temp;

  if (X < Y) {temp = Y; Y = X; X = temp;}

  if (X < TSIL_TOL) return M00000(S);

  if (Y < TSIL_TOL) return M000x0(X,S);

  if (TSIL_FABS(1.0L - Y/X) < TSIL_TOL) return M00xx0(X,S);

  if (TSIL_CABS(S) < TSIL_TOL)
    return (I2(0.0L,X,Y,X) - I2(0.0L,0.0L,X,X) - I2(0.0L,0.0L,Y,X))/(X*Y);

  YoX = Y/X;
  log1mYoX = TSIL_LOG(1.0L - YoX);
  XoY = X/Y;
  logXoY = TSIL_LOG(XoY);
  S = AddIeps(S);
  onemSoY = 1.0L - S/Y;

  if (TSIL_CABS(onemSoY) < 10.0*TSIL_TOL)
    return ((-3.0L*Trilog(Y/(Y-X)) +2.0L*(logXoY+log1mYoX)*Dilog(YoX)
	     -0.5L*log1mYoX*(-log1mYoX*log1mYoX + logXoY*logXoY 
             + 6.0L*Zeta2))/Y);

  onemSoX = 1.0L - S/X;

  if (TSIL_CABS(onemSoX) < 10.0*TSIL_TOL)
    return ((-3.0L*Trilog(1.0L - YoX) +2.0L*(-log1mYoX+I*PI)*Dilog(YoX)
	     +(1.5L*logXoY*log1mYoX - I*PI*logXoY -2.0L*Zeta2)*log1mYoX
	     +2.0L*I*PI*Zeta2)/X);

  log1mSoX = TSIL_CLOG(onemSoX);
  log1mSoY = TSIL_CLOG(onemSoY);
  SmYoSmX = (S-Y)/(S-X);

  return (
  ( log1mSoX*log1mSoX*(-log1mSoX + 2.0L*log1mSoY + log1mYoX + 0.5L*logXoY)  
  + log1mSoY*log1mSoY*(log1mYoX - log1mSoX) 
  + log1mSoX*( (log1mYoX - log1mSoY - 0.5L*logXoY)*logXoY 
      - 2.0L*log1mSoY*log1mYoX - 3.0L*Zeta2 - log1mSoY*TSIL_CLOG(S/X))  
  - log1mSoY*(log1mYoX*logXoY + Zeta2) 
  + (log1mSoY -3.0L*log1mSoX)*Dilog(onemSoX)  
  + (log1mSoX -3.0L*log1mSoY)*Dilog(onemSoY) 
  + (log1mSoY -log1mSoX - 2.0L*logXoY)*Dilog(SmYoSmX)  
  + (log1mSoY -log1mSoX + 2.0L*logXoY)*Dilog(YoX) 
  + 3.0L*((log1mSoY-log1mSoX)*Dilog(XoY*SmYoSmX) 
  + Trilog(onemSoX) + Trilog(onemSoY) - Trilog(SmYoSmX) 
  - Trilog(XoY*SmYoSmX) + Trilog(YoX) - Zeta3))/S);
}

/* ******************************************************************* */
/* D.J. Broadhurst, Z. Phys. C47, 115 (1990), eq. (47) [NOT eq. (48)]  */

TSIL_COMPLEX Mxxyy0 (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S)
{
  TSIL_COMPLEX sqDeltaSXY, tXYS, tYXS;
  TSIL_REAL onemYoX, onemYoX2, onemYoX3, onemYoX4, temp;

  if (X < Y) {temp = Y; Y = X; X = temp;}

  if (X < TSIL_TOL) return M00000(S);

  if (Y < TSIL_TOL) return M00xx0(X,S);
 
  if (TSIL_CABS(S) < TSIL_TOL) {
    onemYoX = 1.0L - Y/X;
    if (TSIL_FABS(onemYoX) > .01)
      return (I2(0.,X,X,X) + I2(0.,Y,Y,X) - 2.0L*I2(0.,X,Y,X))/((X-Y)*(X-Y));
    else {
      onemYoX2 = onemYoX*onemYoX;
      onemYoX3 = onemYoX*onemYoX2;
      onemYoX4 = onemYoX*onemYoX3;

      return ( (1.0L + onemYoX/2.0L + (11.0L/36.0L)*onemYoX2 
		+ (5.0L/24.0L)*onemYoX3 + (137.0L/900.0L)*onemYoX4 
		+ (7.0L/60.0L)*onemYoX4*onemYoX 
		+ (363.0L/3920.0L)*onemYoX4*onemYoX2 
		+ (761.0L/10080.0L)*onemYoX4*onemYoX3 
		+ (7129.0L/113400.0L)*onemYoX4*onemYoX4)/X); 
    }
  }

  if (TSIL_CABS(1.0L - (X+Y+2.0L*TSIL_SQRT(X*Y))/S) < 10.0*TSIL_TOL) {
    TSIL_Warn("Mxxyy0", "M(x,x,y,y,0) is undefined at s = (sqrt(x) + sqrt(y))^2.");
    return TSIL_Infinity;
  }

  S = AddIeps(S);
  sqDeltaSXY = TSIL_CSQRT(Delta(S, X, Y));
  tXYS = ((X+Y-S + sqDeltaSXY)/(2.0L*X));
  tYXS = ((X+Y-S + sqDeltaSXY)/(2.0L*Y));
  return ((F3plusBroadhurst(tXYS) + F3plusBroadhurst(tYXS) 
           - 4.0L*F3plusBroadhurst(TSIL_SQRT(X/Y)*tXYS) 
	   - 4.0L*F3minusBroadhurst(TSIL_SQRT(X/Y)*tXYS) - 6.0L*Zeta3)/S);
}

/* ******************************************************************* */
/* D.J. Broadhurst, Z. Phys. C47, 115 (1990), eq. (33)                 */

TSIL_COMPLEX M0xx0xAtx (TSIL_REAL X)
{
  if (X < TSIL_TOL) {
    TSIL_Warn("M0xx0xAtx", "M(0,x,x,0,x) is undefined for s=x=0.");
    return TSIL_Infinity;
  }

  return (PI*PI*TSIL_LOG(2.0L) - 1.5L*Zeta3)/X;
}

/* ******************************************************************* */
/* Source is unpublished SPM notes (is there a published source?)      */

TSIL_COMPLEX M0xy0yAtx (TSIL_REAL X, TSIL_REAL Y)
{
  TSIL_COMPLEX r, logr, rm1, sqrtr, logrm1;

  if (Y < TSIL_TOL) return M000x0 (X, X);
  if (TSIL_CABS(1.0L - X/Y) < 10.0*TSIL_TOL) return M0xx0xAtx (X);
  if (X < TSIL_TOL) {
    TSIL_Warn("M0xy0yAtx", "M(0,y,x,0,y) is undefined for s=x=0.");
    return TSIL_Infinity;
  }

  r = -AddIeps(-Y/X);
  logr = TSIL_CLOG(r);
  sqrtr = TSIL_CSQRT(r);
  rm1 = r - 1.0L;
  logrm1 = TSIL_CLOG(rm1);

  return (6.0L*Zeta2*TSIL_CLOG(1.0L + sqrtr) - 2.0L*Zeta2*logr 
    - (logrm1*logrm1*logr)/2.0L + logr*logr*logr/3.0L 
    + (logr - 2.0L*logrm1)*Dilog(-rm1) 
    - 2.0L*Trilog((1.0L - sqrtr)/(1.0L + sqrtr)) 
    + 2.0L*Trilog((sqrtr - 1.0L)/(sqrtr + 1.0L)) 
    + 2.0L*Trilog(-rm1) - 2.0L*Trilog(rm1/r) - 1.5L*Zeta3)/X;
}

/* ******************************************************************* */
/* Source is unpublished SPM notes (is there a published source?)      */

TSIL_COMPLEX M0yy0xAtx (TSIL_REAL X, TSIL_REAL Y)
{
  TSIL_COMPLEX r, logr, rm1, rp1, logrm1, logrp1, log2;

  if (Y < TSIL_TOL) return M0000x (X, X);
  if (TSIL_CABS(1.0L - X/Y) < 10.0*TSIL_TOL) return M0xx0xAtx (X);
  if (X < TSIL_TOL) return M0xx00 (Y, 0.0L);

  r = -AddIeps(-Y/X);
  rp1 = r + 1.0L;
  rm1 = r - 1.0L;
  logr = TSIL_CLOG(r);
  logrm1 = TSIL_CLOG(rm1);
  logrp1 = TSIL_CLOG(rp1);
  log2 = 0.693147180559945309417232L;

  return (-4.0L*logr*logr*logr/3.0L + log2*logrp1*logrp1 
         - 2.0L*logrp1*logrp1*logrp1/3.0L 
         - logr*logr*(6.0L*log2 - 3.0L*logrm1 + logrp1)/2.0L 
         + 6.0L*logrp1*Zeta2 - logr*(logrm1*logrm1 - 2.0L*log2*logrp1 
           - 2.0L*logrm1*logrp1 + 6.0L*Zeta2) + Zeta3/2.0L 
         - (log2 + logr/2.0L)*Dilog(1.0L/(r*r)) 
         + 2.0L*(logrm1 - logrp1)*(Dilog(-rm1/2.0L) - Dilog(rm1/(2.0L*r))) 
         + (4.0L*logrm1 - 2.0L*(log2 + logr + logrp1))*Dilog(rm1/r) 
         + (2.0L*log2 + 2.0L*logr - 2.0L*logrp1)*Dilog(r/rp1) 
         - Trilog(1.0L/(r*r))/4.0L - 2.0L*Trilog(rm1/r) 
         - 2.0L*Trilog(r/rp1) - 2*Trilog(-rm1/rp1) 
         + 2.0L*Trilog(rm1/rp1))/X;
}


/* ******************************************************************* */
/* R. Scharf and J.B. Tausk, Nucl. Phys. B412, 523 (1994). eq. (106)   */

TSIL_COMPLEX M0x0y0 (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S)
{
  TSIL_COMPLEX soy, xoy, sox, sqDeltasxy, r1, r2, onemr1, onemr2;
  TSIL_COMPLEX onem1or1, onem1or2, onemsoy, onemsox, onemxoy;
  TSIL_COMPLEX xmsoy, logmsoy, log1msoy, logxoy;
  TSIL_COMPLEX log1msox, logr1, logr2, result;

  if (X < TSIL_TOL) return M000x0(Y,S);
  if (Y < TSIL_TOL) return M000x0(X,S);
  if (TSIL_FABS(1.0L -X/Y) < TSIL_TOL) return M0x0x0 (X, S);
  if (TSIL_CABS(S) < 10.0L*TSIL_TOL) {
    TSIL_Warn("M0x0y0", "M(0,x,0,y,0) is undefined at s = 0.");
    return TSIL_Infinity;
  }

  S = AddIeps(S);

  soy = S/Y; /* Scharf and Tausk's x */
  xoy = X/Y; /* Scharf and Tausk's y */
  sox = S/X;
  sqDeltasxy = TSIL_CSQRT(S*S + X*X + Y*Y - 2.0L*(X*Y + X*S + Y*S));
  r1 = (X+Y-S + sqDeltasxy)/(2.0L*Y);
  r2 = (X+Y-S - sqDeltasxy)/(2.0L*Y);

  onemr1 = 1 - r1;
  onemr2 = 1 - r2;
  onem1or1 = 1 - 1/r1;
  onem1or2 = 1 - 1/r2;
  onemsoy = 1 - soy;
  onemsox = 1 - sox;
  onemxoy = 1 - xoy;
  xmsoy = xoy - soy;
  logmsoy = TSIL_CLOG(-soy);
  log1msoy = TSIL_CLOG(onemsoy);
  log1msox = TSIL_CLOG(onemsox);
  logxoy = TSIL_CLOG(xoy);
  logr1 = TSIL_CLOG(r1);
  logr2 = TSIL_CLOG(r2);

  result = -((logmsoy - logxoy/2.0L)*logr1*logr2
	     + 3.0L*(Trilog(onemr1) + Trilog(onem1or1)
		     + Trilog(onemr2) + Trilog(onem1or2))
	     + 1.5L*logxoy*(Dilog(onemr1) - Dilog(onem1or1) 
			    +Dilog(onemr2) - Dilog(onem1or2))
	     + 0.5L*logxoy*(Dilog(-onemxoy/xoy) - Dilog(onemxoy))
	     -logr1*(Dilog(onemr1) - Dilog(onem1or1))
	     -logr2*(Dilog(onemr2) - Dilog(onem1or2))
	     -0.5L*logxoy*(Dilog(onemr1/onem1or2) - Dilog(onem1or2/onemr1) 
			   + Dilog(onemr2/onem1or1) - Dilog(onem1or1/onemr2))
	     -0.25L*logxoy*logxoy*
             (TSIL_CLOG(onemxoy/onemr2) + TSIL_CLOG(onemxoy/(r1*onemr2))
	      +TSIL_CLOG(onemxoy/onemr1) + TSIL_CLOG(onemxoy/(r2*onemr1)))
	     +(log1msoy*log1msoy*log1msoy + log1msox*log1msox*log1msox)/3.0L
	     + Trilog(soy) + Trilog(sox)
	     + 2.0L*Trilog(-soy/onemsoy) + 2.0L*Trilog(-sox/onemsox)
	     - 2.0L*log1msoy*Dilog(soy) - 2.0L*log1msox*Dilog(sox)
	     -logxoy*(log1msoy*log1msoy - log1msox*log1msox)
	     + (log1msoy + log1msox)*(logr1*logr1 + logr2*logr2)
	     - Trilog(r2*onem1or1) - 2.0L*Trilog(1.0L - r1/onemsoy) 
	     + 2.0L*log1msoy*Dilog(r2*onem1or1)
	     - EtaBranch(r2, onemsoy)*
             (0.5L*TSIL_CLOG(r2*onem1or1)*TSIL_CLOG(r2*onem1or1) 
	      - TSIL_CLOG(1 - r1/onemsoy)*TSIL_CLOG(1 - r1/onemsoy) 
	      - 2.0L*log1msoy*TSIL_CLOG(r2*onem1or1))
	     - Trilog(onemr1/r2) -2.0L*Trilog(1.0L - r2/xmsoy) 
	     + 2.0L*log1msox*Dilog(onemr1/r2)
	     - EtaBranch(r1, xmsoy)*
             (0.5L*TSIL_CLOG(onemr1/r2)*TSIL_CLOG(onemr1/r2) 
	      - TSIL_CLOG(1-r2/xmsoy)*TSIL_CLOG(1-r2/xmsoy) 
	      - 2.0L*log1msox*TSIL_CLOG(onemr1/r2))
	     - Trilog(r1*onem1or2) - 2.0L*Trilog(1.0L - r2/onemsoy) 
	     + 2.0L*log1msoy*Dilog(r1*onem1or2)
	     - EtaBranch(r1, onemsoy)*
             (0.5L*TSIL_CLOG(r1*onem1or2)*TSIL_CLOG(r1*onem1or2) 
	      - TSIL_CLOG(1.0L - r2/onemsoy)*TSIL_CLOG(1.0L - r2/onemsoy) 
	      - 2.0L*log1msoy*TSIL_CLOG(r1*onem1or2))
	     - Trilog(onemr2/r1) -2.0L*Trilog(1.0L - r1/xmsoy) 
	     + 2.0L*log1msox*Dilog(onemr2/r1)
	     - EtaBranch(r2, xmsoy)*
             (0.5L*TSIL_CLOG(onemr2/r1)*TSIL_CLOG(onemr2/r1) 
	      - TSIL_CLOG(1.0L-r1/xmsoy)*TSIL_CLOG(1.0L-r1/xmsoy) 
	      - 2.0L*log1msox*TSIL_CLOG(onemr2/r1)))/S;

  return result;
}

/* ******************************************************************* */
/* R. Scharf and J.B. Tausk, Nucl. Phys. B412, 523 (1994). eq. (105)   */

TSIL_COMPLEX M0x0x0 (TSIL_REAL X, TSIL_COMPLEX S)
{
  TSIL_COMPLEX rat, rat2, r1, r2, onemrat, onemr1; 
  TSIL_COMPLEX onemr2, logmrat, logonemrat, onemr1mr2;

  if (X < TSIL_TOL) return M00000 (S);

  S = AddIeps(S);
  rat = S/X;

  if (TSIL_CABS(rat) < 10.0*TSIL_TOL) {
    TSIL_Warn("M0x0x0", "M(0,x,0,x,0) is undefined at s = 0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(rat) < 0.001L) {
    logmrat = TSIL_CLOG(-rat);
    rat2 = rat*rat;
 
    return (2.0L - logmrat - (1.0L + logmrat)*rat/12.0L 
	    - (23.0L + 6.0L*logmrat)*rat2/540.0L 
	    - (17.0L + 2.0L*logmrat)*rat2*rat/1120.0L 
	    - (367.0L + 20.0L*logmrat)*rat2*rat2/63000.0L 
	    - (2531.0L + 60.0L*logmrat)*rat2*rat2*rat/997920.0L)/X;
  }

  if (TSIL_CABS(1.0L - rat/4.0L) < TSIL_TOL)
    return (-2.41221075861979707519227L + 3.33855021948523315348115L*I)/X;

  if (TSIL_CABS(1.0L - rat) < TSIL_TOL)
    return (1.847263653039733073915529L + 3.445141853366646686164035L*I)/X;

  onemrat = 1.0L -rat;
  logonemrat = TSIL_CLOG(onemrat);

  r1 = 1.0L - (S - TSIL_CSQRT(S*(S-4.0L*X)))/(2.0L*X);
  r2 = 1.0L - (S + TSIL_CSQRT(S*(S-4.0L*X)))/(2.0L*X);
  onemr1 = 1.0L -r1;
  onemr2 = 1.0L -r2;
  onemr1mr2 = 1.0L -r1-r2;

  return -((2.0L*logonemrat*logonemrat*logonemrat)/3.0L
  + EtaBranch(onemrat, r1)*
    (2.0L*TSIL_CLOG(onemr1/onemr1mr2)*TSIL_CLOG(onemr1/onemr1mr2)
    - TSIL_CLOG(onemr1*r1)*TSIL_CLOG(onemr1*r1))
  + EtaBranch(onemrat, r2)*
    (2.0L*TSIL_CLOG(onemr2/onemr1mr2)*TSIL_CLOG(onemr2/onemr1mr2)
    - TSIL_CLOG(onemr2*r2)*TSIL_CLOG(onemr2*r2))
  + TSIL_CLOG(r2)*TSIL_CLOG(r2)*(-2*logonemrat - TSIL_CLOG(-rat))
  + 2.0L*(Dilog(onemr1) - Dilog(onemr2))*TSIL_CLOG(r2)
  + 6.0L*Trilog(onemr1) + 6.0L*Trilog(onemr2)
  - 4.0L*Trilog(onemr1/onemr1mr2) - 4.0L*Trilog(onemr2/onemr1mr2)
  - 2.0L*Trilog(onemr1*r1) - 2.0L*Trilog(onemr2*r2)
  + 2.0L*Trilog(rat) + 4.0L*Trilog(-rat/onemrat))/S;
}

/* ******************************************************************* */
/* J. Fleischer, A.V. Kotikov and O.L. Veretin,                        */
/* Nucl. Phys. B547, 343 (1999) [hep-ph/9808242] I_15, eq. (77)        */

TSIL_COMPLEX M000xx (TSIL_REAL X, TSIL_COMPLEX S)
{
  TSIL_COMPLEX z, y, logy, r, logr, r2, r3, r4, r5;

  if (X < TSIL_TOL) return M00000(S);

  S = AddIeps(S);
  z = S/X; 
             
  if (TSIL_CABS(z) < 10.0L*TSIL_TOL) {
    TSIL_Warn("M000xx", "M(0,0,0,x,x) is undefined for s=0.");
    return TSIL_Infinity;
  }

  r = 1.0L - z;
  logr = TSIL_CLOG(r);      

  if (TSIL_CABS(r) < 0.002L) {
    r2 = r*r;
    r3 = r2*r;
    r4 = r3*r;
    r5 = r4*r;
             
    return ((2.77089547955959961087329537L + 5.1677127800499700292460517628L*I)
     + r*(-2.02438012507831050117489046L + PI*I*(logr - 1.0L))
     + r2*(-0.3087402214806007879389257L + PI*I*(logr/2.0L - 1.0L/4.0L))
     + r3*(-0.1328379607518328669262794L + PI*I*(logr/3.0L - 1.0L/9.0L))
     + r4*(-0.0708048756412441030731011L + PI*I*(logr/4.0L - 1.0L/16.0L))
     + r5*(-0.0442881314475634082095990L + PI*I*(logr/5.0L - 1.0L/25.0L)))/S;
  }

  y = (TSIL_CSQRT(z-4.0L) - TSIL_CSQRT(z))/(TSIL_CSQRT(z-4.0L) + TSIL_CSQRT(z));
  logy = TSIL_CLOG(y);

  return (2.0L*Trilog(z) - TSIL_CLOG(-z)*Dilog(z) + Zeta2*logr
         + logy*logy*logy/6.0L - 0.5L*logy*logy*
            (8.0L*TSIL_CLOG(1.0L - y) - 3.0L*TSIL_CLOG(1.0L - y + y*y))
         - 6.0L*Zeta3 - Trilog(-y*y*y)/3.0L
         + 3.0L*Trilog(-y) + 8.0L*Trilog(y) +
         + logy*(Dilog(-y*y*y) - 3.0L*Dilog(-y) - 8.0L*Dilog(y)))/S;
}

/* ******************************************************************* */     
/* J. Fleischer, A.V. Kotikov and O.L. Veretin,                        */
/* Nucl. Phys. B547, 343 (1999) [hep-ph/9808242] I_14, eq. (76)        */

TSIL_COMPLEX M0xx00 (TSIL_REAL X, TSIL_COMPLEX S)
{
  TSIL_COMPLEX z, logmz, y, logy, r, logr, logr2, r2, r3, r4, r5;

  if (X < TSIL_TOL) return M00000(S);
             
  S = AddIeps(S);
  z = S/X; 

  if (TSIL_CABS(z) < TSIL_TOL) return 2.0L*Zeta2/X;

  logmz = TSIL_CLOG(-z);

  if (TSIL_CABS(z) < 0.000001L) 
     return (2.0L*Zeta2 + z*(0.25L + Zeta2 - 0.5L*logmz)
            + z*z*(-13.0L/36.0L + 2.0L*Zeta2/3.0L - logmz/3.0L))/X; 
           
  r = 1.0L - z;
  logr = TSIL_CLOG(r);      
  r2 = r*r;

  if ( TSIL_CABS(r) < 0.002L ) {
    r3 = r2*r;
    r4 = r3*r;
    r5 = r4*r;
    logr2 = logr*logr;                   

    return (6.37706618903838246707250487L + I*PI*PI*PI/6.0L  
      + r*(-9.11071228462919048942954040L - 2.0L*PI*I + 2.0L*PI*I*logr) 
      + r2*(4.64845646019748631941773107L + 0.5L*PI*I 
         + (1.0L - PI*I)*logr - logr2) 
      + r3*(-3.2449486806021261963444960L - 2.0L*PI*I/9.0L  
         + 2.0L*PI*I*logr/3.0L) 
      + r4*(2.42467543171744466412626622L + 0.125L*PI*I  
         + (0.25L - 0.5L*PI*I)*logr - 0.5*logr2) 
      + r5*(-1.9644518835048194797987766L - 0.08L*PI*I + 0.4L*PI*I*logr))/S;
  }

  y = (TSIL_CSQRT(z-4.0L) - TSIL_CSQRT(z))/(TSIL_CSQRT(z-4.0L) + TSIL_CSQRT(z));
  logy = TSIL_CLOG(y);

  return (2.0L*S12FKV(1.0L/r) +2.0L*S12FKV(-1.0L/r) -S12FKV(1.0L/r2) 
           + TSIL_CLOG(2.0L - z)*(logr*logr - 2.0L*logr*logmz 
             - 2.0L*Dilog(z)) - 2.0L*logr*logr*logr/3.0L 
           + logr*logr*logmz - 2.0L*Zeta2*logr + logy*logy*logy/3.0L
           + logy*logy*(2.0L*TSIL_CLOG(1.0L+y*y) - 3.0L*TSIL_CLOG(1.0L-y+y*y)) 
           - 6.0L*Zeta3 -Trilog(-y*y) + 2.0L*Trilog(-y*y*y)/3.0L 
           - 6.0L*Trilog(-y) + 2.0L*logy*(Dilog(-y*y) 
           - Dilog(-y*y*y) + 3.0L*Dilog(-y)))/S;
}

/* ******************************************************************* */
/* J. Fleischer, A.V. Kotikov and O.L. Veretin,                        */
/* Nucl. Phys. B547, 343 (1999) [hep-ph/9808242] I_123, eq. (78)       */

TSIL_COMPLEX M0xxx0 (TSIL_REAL X, TSIL_COMPLEX S)
{
  TSIL_COMPLEX z, y, logy, r, logr, r2, r3, r4, r5;

  if (X < TSIL_TOL) return M00000(S);
             
  S = AddIeps(S);
  z = S/X; 

  if (TSIL_CABS(z) < 0.002L) return (Zeta2 + z*(0.5L*Zeta2 - 0.25L) 
            + z*z*(Zeta2/3.0L - 7.0L/24.0L)
            + z*z*z*(Zeta2/4.0L -131.0L/480.0L)
            + z*z*z*z*(Zeta2/5.0L -2459.0L/10080.0L)
            + z*z*z*z*z*(Zeta2/6.0L -9367.0L/43200.0L)
            + z*z*z*z*z*z*(Zeta2/7.0L -107417.0L/554400.0L))/X; 
           
  r = 1.0L - z;
  logr = TSIL_CLOG(r);      
  r2 = r*r;

  if( TSIL_CABS(r) < 0.002 )
  {
    r3 = r2*r;
    r4 = r3*r;
    r5 = r4*r;             
      
    return (3.04932055619932735931526231457262909045L 
             + r*(-2.45564510912370625588737537084817293002L 
               + 1.81379936423421785059407825764215573229L*logr) 
             + r2*(-0.35697429081491473421554951886766953296L 
               + 0.30229989403903630843234637627369262205L*logr) 
             + r3*(-0.08747177297904619871428937340651466152L 
               + 0.20153326269269087228823091751579508137L*logr) 
             + r4*(-0.05110761517969126259235005454346829001L 
               + 0.1175610699040696755014680352175471308L*logr) 
             + r5*(-0.02550680592748196593530901978714506521L 
               + 0.08509182202580281274391972072889125658L*logr))/S;
  }

  y = (TSIL_CSQRT(z-4.0L) - TSIL_CSQRT(z))/(TSIL_CSQRT(z-4.0L) + TSIL_CSQRT(z));
  logy = TSIL_CLOG(y);

  return (-Zeta2*(logr + logy) - 6.0L*Zeta3 
           - 1.5L*TSIL_CLOG(1.0L - y + y*y)*logy*logy 
           + Trilog(-y*y*y) - 9.0L*Trilog(-y) 
           - 2.0L*logy*(Dilog(-y*y*y) - 3.0L*Dilog(-y)))/S;
}

/* ******************************************************************* */
/* J. Fleischer, A.V. Kotikov and O.L. Veretin,                        */
/* Nucl. Phys. B547, 343 (1999) [hep-ph/9808242] I_14, eq. (69)        */

TSIL_COMPLEX S12FKV (TSIL_COMPLEX z)
{
  TSIL_COMPLEX onemz, onemz2, log1mz, log1mz2;

  if (TSIL_CABS(z) < .002L) return z*z*(0.25L + z/6.0L + 11.0L*z*z/96.0L 
     + z*z*z/12.0L + 137.0L*z*z*z*z/2160.0L + z*z*z*z*z/20.0L);    
  
  onemz = 1.0L - z;
  onemz2 = onemz*onemz;
  log1mz = TSIL_CLOG(onemz);
  log1mz2 = log1mz*log1mz;

  if (TSIL_CABS(onemz) < 0.003L)
    return (Zeta3 
     + (-1.0L + log1mz - log1mz2/2)*onemz 
     + (-1.0L/8.0L + log1mz/4.0L - log1mz2/4.0L)*onemz2 
     + (-1.0L/27.0L + log1mz/9.0L - log1mz2/6.0L)*onemz2*onemz 
     + (-1.0L/64.0L + log1mz/16.0L - log1mz2/8.0L)*onemz2*onemz2 
     + (-1.0L/125.0L + log1mz/25.0L - log1mz2/10.0L)*onemz2*onemz2*onemz
     + (-1.0L/216.0L + log1mz/36.0L - log1mz2/12.0L)*onemz2*onemz2*onemz2);

  return 0.5L*log1mz2*TSIL_CLOG(z) + log1mz*Dilog(onemz) - Trilog(onemz) + Zeta3;
}

/* ***************************************************************** */
 
TSIL_COMPLEX sMAtZero (TSIL_REAL x,
		       TSIL_REAL y,
		       TSIL_REAL z,
		       TSIL_REAL u,
		       TSIL_REAL v)
{
  return 0.0L + 0.0L*I;
}
  
/* ***************************************************************** */
            
TSIL_COMPLEX sMprimeAtZero (TSIL_REAL x,
			    TSIL_REAL y, 
			    TSIL_REAL z,
			    TSIL_REAL u,
			    TSIL_REAL v)
{
  TSIL_REAL tmpx, tmpz, qq;
                     
  /* Result should be independent of qq, but need to pick something. */
  qq = x + y + z + u + v;
                     
  if (TSIL_FABS(x-z) > TSIL_FABS(y-u))
    {tmpx = y; tmpz = u; y = x; u = z; x = tmpx; z = tmpz;}
 
  if (TSIL_FABS(z - x) > TSIL_TOL)
    return (I2(x, y, v, qq) - I2(x, u, v, qq)
            - I2(z, y, v, qq) + I2(z, u, v, qq))/((x - z)*(y - u));

  else if (TSIL_FABS(y - u) > TSIL_TOL)
    return (I2p(x, y, v, qq) - I2p(x, u, v, qq))/(y - u);
                          
  else
    return I2pp(x, y, v, qq);
}
