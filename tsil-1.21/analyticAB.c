/* Analytic results for A-type and B-type one-loop functions */

#include "internal.h"

/* **************************************************************** */
/* Wrappers for the user API functions:                             */

TSIL_COMPLEX TSIL_A (TSIL_REAL X, TSIL_REAL QQ)
{
  return A(X, QQ);
}

TSIL_COMPLEX TSIL_Ap (TSIL_REAL X, TSIL_REAL QQ)
{
  return A(X, QQ);
}

TSIL_COMPLEX TSIL_Aeps (TSIL_REAL x, TSIL_REAL qq)
{
  return Aeps (x, qq);
}

TSIL_COMPLEX TSIL_B (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  return B (X, Y, S, QQ);
}

TSIL_COMPLEX TSIL_Bp (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{  
  return Bp (X, Y, S, QQ);
}

TSIL_COMPLEX TSIL_dBds (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{  
  return dBds_analytic (X, Y, S, QQ);
}

TSIL_COMPLEX TSIL_Beps (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  return Beps (X, Y, S, QQ);
}

/* ******************************************************************* */
/* Implements Eq. 2.10 of hep-ph/0307101                               */
/* Modified in v1.2 to accept negative squared mass arg                */

TSIL_COMPLEX A (TSIL_REAL x, TSIL_REAL qq)
{
  if (TSIL_FABS(x) < TSIL_TOL) 
    return 0.0;
  if (x > 0) 
    return (x * (TSIL_LOG(x/qq) - 1.));  
  return (x * (TSIL_LOG(-x/qq) - 1. + I*PI));
}

/* ******************************************************************* */
/* A'(x), added to TSIL API in v1.2                                    */

TSIL_COMPLEX Ap (TSIL_REAL x, TSIL_REAL qq)
{
  if (TSIL_FABS(x) < TSIL_TOL) 
    return 0.0;
  if (x > 0) 
    return (TSIL_LOG(x/qq));  
  return (TSIL_LOG(-x/qq) + I*PI);
}

/* ******************************************************************** */

TSIL_COMPLEX Aeps (TSIL_REAL x, TSIL_REAL qq)
{
  TSIL_COMPLEX lnbarx = Ap (x, qq);
  
  if (TSIL_FABS(x) < TSIL_TOL) return 0.0L;
  else return x * (-1.0L - 0.5L*Zeta2 + lnbarx - 0.5L*lnbarx*lnbarx);
}

/* **************************************************************** */
/* hep-ph/0307101 eq. (6.3)                                         */

TSIL_COMPLEX B (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_REAL    temp;
  TSIL_COMPLEX sqDeltaSXY, lnbarX, lnbarY;

  if (TSIL_FABS (X) < TSIL_FABS (Y)) {temp = Y; Y = X; X = temp;}

  if (TSIL_FABS (X) < TSIL_TOL)
    return B00(S,QQ);

  if (TSIL_FABS (Y) < TSIL_TOL)
    return B0x(X,S,QQ);

  if (TSIL_CABS (S) < TSIL_TOL) {
    if (TSIL_FABS (1.0L - Y/X) > 0.0L)
      return (1.0L + (Y*TSIL_LOG(Y/QQ) - X*TSIL_LOG(X/QQ))/(X-Y));
    else 
      return (-TSIL_LOG (X/QQ));
  }

  S = AddIeps(S);   
  sqDeltaSXY = TSIL_CSQRT(Delta(S, X, Y));
  lnbarX = TSIL_LOG (X/QQ);
  lnbarY = TSIL_LOG (Y/QQ);

  /* Following avoids roundoff error for very negative s. */
  if ((TSIL_CREAL(S) < -10.0L*(X+Y)) && (TSIL_CIMAG(S) < TSIL_TOL)) {   
    return (2.0L - 0.5L * (lnbarX + lnbarY) +
	  (sqDeltaSXY * TSIL_CLOG(0.5L*(X + Y - S + sqDeltaSXY)/Y) +
	   0.5L * (Y - X - sqDeltaSXY) * (lnbarX - lnbarY))/S);
  }

  return (2.0L - 0.5L * (lnbarX + lnbarY) +
	  (-sqDeltaSXY * TSIL_CLOG(0.5L*(X + Y - S - sqDeltaSXY)/X) +
	   0.5L * (Y - X - sqDeltaSXY) * (lnbarX - lnbarY))/S);
}

/* **************************************************************** */
/* hep-ph/0307101 eq. (6.5)                                         */

TSIL_COMPLEX B00 (TSIL_COMPLEX S, TSIL_REAL QQ)
{
  if (TSIL_CABS (S) < TSIL_TOL) {
    TSIL_Warn("B00", "B(0,0) is undefined when s=0."); 
    return TSIL_Infinity;
  }

  S = AddIeps(S);    
  return (2.0L - TSIL_CLOG(-S/QQ));
}

/* **************************************************************** */
/* hep-ph/0307101 eq. (6.4)                                         */

TSIL_COMPLEX B0x (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  if (TSIL_FABS (X) < TSIL_TOL)
    return B00(S,QQ);

  if (TSIL_CABS (S) < TSIL_TOL)
    return (1.0L - TSIL_LOG (X/QQ));

  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL)
    return 2.0L - TSIL_LOG(X/QQ);

  S = AddIeps(S);  
  return 2.0L + ((X - S)*TSIL_CLOG((X - S)/QQ) - X*TSIL_LOG(X/QQ))/S;
}

/* ******************************************************************** */
/* Following is the derivative of B with respect to its first argument. */
/* hep-ph/0307101 eq. (3.1)                                             */

TSIL_COMPLEX Bp (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  if (X < TSIL_TOL) {
    TSIL_Warn("Bp", "B(x',y) is undefined for x=0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(1.0L - S/(X+Y+2.0L*TSIL_SQRT(X*Y))) < TSIL_TOL) {
    TSIL_Warn("Bp", "B(x',y) is undefined at s = (sqrt(x) + sqrt(y))^2.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(S) < TSIL_TOL) {
    if (TSIL_FABS(1.0L - X/Y) < TSIL_TOL)
      return (-0.5L/X);
    else
      return 1.0L/(Y-X) + Y*TSIL_LOG(X/Y)/((Y-X)*(Y-X));
  }

  if (TSIL_CABS(1.0L - (X + Y - 2.0L*TSIL_SQRT(X*Y))/S) < TSIL_TOL) 
    return (1.0L - TSIL_SQRT(Y/X) +0.5L*TSIL_LOG(Y/X))/(X + Y - 2.0L*TSIL_SQRT(X*Y));
  else
    return ((X-Y-S)*B(X,Y,S,QQ) + (X+Y-S)*TSIL_LOG(X/QQ) 
	    -2.0L*A(Y,QQ) + 2.0L*(S-X))/Delta(S,X,Y);
}

/* ******************************************************************** */
/* Following is derivative of B(x,y) with respect to s.                 */
/* Note it is quite distinct from dBds() used in the Runge-Kutta!       */
/* This one is a stand-alone analytic function.                         */

TSIL_COMPLEX dBds_analytic (TSIL_REAL x,
			    TSIL_REAL y,
			    TSIL_COMPLEX s,
			    TSIL_REAL qq)
{
  TSIL_REAL thxy, psxy, sqrtx, sqrty, alphax, alphay, temp;
  TSIL_COMPLEX Btemp, result;

  if (x < y) {temp = y; y = x; x = temp;}

  if (x + TSIL_CABS(s) < TSIL_TOL) {
    TSIL_Warn("dBds_analytic", "dBds(0,0) is undefined when s=0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(s/x) < TSIL_TOL) return BprimeAtZero (x, y, qq);

  thxy = Th2 (x, y);
  psxy = Ps2 (x, y);
  sqrtx = TSIL_SQRT(x);
  sqrty = TSIL_SQRT(y);
  alphax = Alpha(x,qq);
  alphay = Alpha(y,qq);
  
  if (TSIL_CABS(1.0L - s/thxy) < TSIL_TOL) {
    TSIL_Warn("dBds_analytic", "dBds(x,y) is undefined at threshold.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS(1.0L - s/psxy) < TSIL_TOL) {
    return (0.5L * alphax * (1.0L + sqrty/sqrtx) 
	    - 0.5L * alphay * (1.0L + sqrtx/sqrty)
	    + 2.0L * (sqrty - sqrtx) )/TSIL_POW(sqrtx - sqrty,3);
  }

  Btemp = B (x, y, s, qq);

  result = (0.5L/s) * (
           (thxy * Btemp - s + (sqrtx + sqrty)*(alphax + alphay))/(s - thxy)
         + (psxy * Btemp - s + (sqrtx - sqrty)*(alphax - alphay))/(s - psxy));
 
  return result;
}

/* ******************************************************************** */
/* R. Scharf and J.B. Tausk, Nucl. Phys. B412, 523 (1994). eq. (83)     */

TSIL_COMPLEX BepsAtZero (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq)
{
  TSIL_COMPLEX lnbarx, lnbary;
  TSIL_REAL temp;

  if (x < y) {temp = x; x = y; y = temp;}

  if (x < TSIL_TOL) {
    TSIL_Warn("BepsAtZero", "Beps(0,0) is undefined at s = 0.");
    return TSIL_Infinity;
  }

  lnbarx = TSIL_CLOG(x/qq);

  if (y < TSIL_TOL) return 1.0L + Zeta2/2.0L - lnbarx + lnbarx*lnbarx/2.0L;

  if (TSIL_FABS(x-y)/(x+y) < TSIL_TOL) return (Zeta2 + lnbarx*lnbarx)/2.0L;

  lnbary = TSIL_CLOG(y/qq);

  return 1.0L + Zeta2/2.0L + (x*lnbarx*(lnbarx/2.0L - 1.0L) -
			      y*lnbary*(lnbary/2.0L - 1.0L))/(x-y);
}

/* ******************************************************************** */
/* R. Scharf and J.B. Tausk, Nucl. Phys. B412, 523 (1994). eq. (83)     */

TSIL_COMPLEX Beps0x (TSIL_REAL x, TSIL_COMPLEX s, TSIL_REAL qq)
{
  TSIL_COMPLEX sqrtx, lnbarx, log1msox, lnbarms;

  if (TSIL_CABS(s) < TSIL_TOL) return BepsAtZero (0, x, qq);

  if (x < TSIL_TOL) {
    lnbarms = TSIL_CLOG(-AddIeps(s)/qq);
    return 4.0L - Zeta2/2.0L - 2.0L*lnbarms + lnbarms*lnbarms/2.0L;
  }

  sqrtx = TSIL_CSQRT(x);
  lnbarx = TSIL_CLOG(x/qq);

  if (TSIL_CABS(1 - s/x) < 10.0*TSIL_TOL)
    return 4.0L + Zeta2/2.0L + lnbarx*(lnbarx - 4.0L)/2.0L;

  log1msox = TSIL_CLOG(1.0L - AddIeps(s)/x);

  return 4.0L + Zeta2/2.0L - 2.0L*lnbarx + lnbarx*lnbarx/2.0L 
         + (1.0L - x/s)*(lnbarx*log1msox - 2.0L*log1msox 
         + log1msox*log1msox/2.0L -Dilog(s/(s-x)));
}

/* ******************************************************************** */
/* R. Scharf and J.B. Tausk, Nucl. Phys. B412, 523 (1994). eq. (83)     */

TSIL_COMPLEX Beps (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX sqrtdeltasxy, t1, t2, t3, t4, logt1, logt2, logt3, logt4;
  TSIL_COMPLEX log1mt1, log1mt2, sqrtx, sqrty, lnbarx, lnbary;
  TSIL_REAL temp;

  if (X < Y) {temp = X; X = Y; Y = temp;}

  if (TSIL_CABS(S) < TSIL_TOL) return BepsAtZero (X, Y, QQ);

  if (Y < TSIL_TOL) return Beps0x (X, S, QQ);
  
  sqrtx = TSIL_CSQRT(X);
  sqrty = TSIL_CSQRT(Y);
  lnbarx = TSIL_CLOG(X/QQ);
  lnbary = TSIL_CLOG(Y/QQ);

  if (TSIL_CABS(S - (sqrtx + sqrty)*(sqrtx + sqrty))/(X+Y) < 10.0*TSIL_TOL) 
    return 4.0L + Zeta2/2.0L 
           + (sqrtx*lnbarx*(lnbarx - 4.0L) + 
	      sqrty*lnbary*(lnbary - 4.0L))/(2.0L*(sqrtx + sqrty)); 

  if (TSIL_CABS(S - (sqrtx - sqrty)*(sqrtx - sqrty))/(X+Y) < 10.0*TSIL_TOL) 
    return 4.0L + Zeta2/2.0L 
           + (sqrtx*lnbarx*(lnbarx - 4.0L) - 
	      sqrty*lnbary*(lnbary - 4.0L))/(2.0L*(sqrtx - sqrty)); 

  S = AddIeps(S);

  sqrtdeltasxy = TSIL_CSQRT(Delta(S,X,Y));

  t1 = ( S - X + Y + sqrtdeltasxy)/(2.0L*sqrtdeltasxy);
  t2 = (-S - X + Y + sqrtdeltasxy)/(2.0L*sqrtdeltasxy);
  t3 = (-S + X + Y + sqrtdeltasxy)/(2.0L*X);
  t4 = (-S + X + Y - sqrtdeltasxy)/(2.0L*X);

/*   printf("\ns = ");TSIL_cprintf(S);printf("\n"); */
/*   printf("x = %Lf\n", X); */
/*   printf("y = %Lf\n", Y); */
/*   printf("sqrtdeltasxy = ");TSIL_cprintf(sqrtdeltasxy);printf("\n"); */

/*   printf("\nt1    = ");TSIL_cprintf(t1);printf("\n"); */
/*   printf("t2    = ");TSIL_cprintf(t2);printf("\n"); */
/*   printf("t3    = ");TSIL_cprintf(t3);printf("\n"); */
/*   printf("t4    = ");TSIL_cprintf(t4);printf("\n"); */
  
  logt1 = TSIL_CLOG(t1);
  logt2 = TSIL_CLOG(t2);
  logt3 = TSIL_CLOG(t3);
  logt4 = TSIL_CLOG(t4);

/*   printf("\nlogt1 = ");TSIL_cprintf(logt1);printf("\n"); */
/*   printf("logt2 = ");TSIL_cprintf(logt2);printf("\n"); */
/*   printf("logt3 = ");TSIL_cprintf(logt3);printf("\n"); */
/*   printf("logt4 = ");TSIL_cprintf(logt4);printf("\n"); */

  log1mt1 = TSIL_CLOG(1.0L - t1);
  log1mt2 = TSIL_CLOG(1.0L - t2);

  return 4.0L + Zeta2/2.0L + (lnbarx*lnbarx + lnbary*lnbary)/4.0L 
         - lnbarx - lnbary + (sqrtdeltasxy*(Dilog(t2) - Dilog(t1) 
              + (logt3 - logt4)*(1.0L - lnbarx/4.0L - lnbary/4.0L) 
              + (log1mt1 + log1mt2)*(logt2 - logt1)/2.0L) 
         + (X - Y)*(lnbary - lnbarx)*(1.0L - lnbarx/4.0L - lnbary/4.0L))/S;
}

/* **************************************************************** */
  
TSIL_COMPLEX BAtZero (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq)
{
  if (TSIL_FABS(x - y) > TSIL_TOL)
    return (A(y,qq) - A(x,qq))/(x - y);
  else  
    return -A(x,qq)/x - 1.0L;
}
    
/* **************************************************************** */
    
TSIL_COMPLEX BprimeAtZero (TSIL_REAL x, TSIL_REAL y, TSIL_REAL qq)
{
  TSIL_REAL xoy, onemxoy, onemxoy2, onemxoy3, onemxoy4, temp;
              
  if (x < y) {temp = y; y = x; x = temp;}

  xoy = x/y;
  onemxoy = 1.0L - xoy;
 
  if (TSIL_FABS(onemxoy) > 0.005) 
    return (x*x - 2.L*A(x,qq)*y - y*y + 2.L*x*A(y,qq))/(2.L*TSIL_POW(x - y,3));
  else {
    onemxoy2 = onemxoy * onemxoy;
    onemxoy3 = onemxoy2 * onemxoy;
    onemxoy4 = onemxoy3 * onemxoy;

    return (1.0L - onemxoy/2.0L - onemxoy2/5.0L - onemxoy3/10.0L
              - 2.0L*onemxoy4/35.0L - onemxoy2*onemxoy3/28.0L
              - onemxoy3*onemxoy3/42.0L - onemxoy4*onemxoy3/60.0L
              - 2.0L*onemxoy4*onemxoy4/165.0L
              - onemxoy3*onemxoy3*onemxoy3/110.0L
              - onemxoy3*onemxoy3*onemxoy4/143.0L
              - onemxoy3*onemxoy4*onemxoy4/182.0L
              - 2.0L*onemxoy4*onemxoy4*onemxoy4/455.0L)/(6.0L * x);
  }
}
