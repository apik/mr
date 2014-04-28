/* Analytic results for V-type functions. */

#include "internal.h"

/* ******************************************************************* */
/* A wrapper for the user API:                                         */

int TSIL_Vanalytic (TSIL_REAL X,
		    TSIL_REAL Y,
		    TSIL_REAL Z,
		    TSIL_REAL U, 
		    TSIL_COMPLEX S,
		    TSIL_REAL QQ,
		    TSIL_COMPLEX *result)
{
  return Vanalytic (X, Y, Z, U, S, QQ, result);
}

/* ******************************************************************* */
/* Defined in hep-ph/0307101 eq. (2.22)                                */

int Vanalytic (TSIL_REAL X,
	       TSIL_REAL Y,
	       TSIL_REAL Z,
	       TSIL_REAL U, 
               TSIL_COMPLEX S,
	       TSIL_REAL QQ,
	       TSIL_COMPLEX *result)
{
  int success = 1;

  if (Y < TSIL_TOL) {
    TSIL_Warn("Vanalytic", "V(x,y,z,u) is undefined for y = 0.");
    *result = TSIL_Infinity;
  }
  else if (TSIL_CABS(S) < TSIL_TOL) 
    {
      if (TSIL_FABS(1.0L - Y/X) > TSIL_TOL)
	*result = ( (I2(X, Z, U, QQ) - I2(Y, Z, U, QQ) + 
		     (Y-X)*I2p(Y, Z, U, QQ))/((X-Y)*(X-Y)) );
      else *result = 0.5L*I2p2 (X, Z, U, QQ);
    }
  else if ((X < TSIL_TOL) && (TSIL_FABS(Delta(Y,Z,U))/(Y*Y+Z*Z+U*U) > TSIL_TOL)) 
    {
      *result = V0xyz(Y, Z, U, S, QQ);
    }
  else if (U + TSIL_FABS(Y-Z) < TSIL_TOL) *result = Vxy0y(X,Y,S,QQ);
  else if (Z + TSIL_FABS(Y-U) < TSIL_TOL) *result = Vxy0y(X,Y,S,QQ);
  else success = 0;

  return success;
}

/* ******************************************************************* */
/* Special case of hep-ph/0307101 eq. (3.22)                           */

TSIL_COMPLEX V0x00 (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_REAL lnbarX;  

  if (X < TSIL_TOL) {
    TSIL_Warn ("V0x00", "V(0,x,0,0) is undefined for x=0.");
    return TSIL_Infinity;
  }

  lnbarX = TSIL_LOG(X/QQ);

  if (TSIL_CABS(S) < TSIL_TOL)
    return (2.0L - lnbarX)/X;

  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL) {
    TSIL_Warn("V0x00", "V(0,x,0,0) is undefined for s=x.");
    return TSIL_Infinity;
  }

  return ((X*U0x00(X,S,QQ) -2.0L*S000(S,QQ) +I200x(X,QQ) + 1.25L*S -X)/(S-X) 
	  + B0x(X,S,QQ))/X;
}

/* ******************************************************************* */
/* hep-ph/0307101 eq. (6.21)                                           */

TSIL_COMPLEX Vxy0y (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_REAL lnbarY;

  if (Y < TSIL_TOL) {
    TSIL_Warn("Vxy0y", "V(x,y,0,y) is undefined for y=0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS (1.0L - S/(X + Y + 2.0L*TSIL_SQRT(X*Y))) < 10.0L*TSIL_TOL) {
    TSIL_Warn("Vxy0y", "V(x,y,0,y) is undefined at s = (sqrt(x)+sqrt(y))^2.");
    return TSIL_Infinity;
  }

  lnbarY = TSIL_LOG(Y/QQ);

  return (Tbar0xy(X,Y,S,QQ) - Tx0y(Y,X,S,QQ) -lnbarY*B(X,Y,S,QQ))/(2.0L*Y)
	 + (lnbarY - 2.0L)*Bp(Y,X,S,QQ);
}

/* ******************************************************************* */
/* Special case of hep-ph/0307101 eq. (6.21)                           */

TSIL_COMPLEX V0x0x (TSIL_REAL X, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_COMPLEX log1mSoX;
  TSIL_REAL    lnbarX;

  if (TSIL_FABS(X) < TSIL_TOL) {
    TSIL_Warn("V0x0x", "V(0,x,0,x) is undefined for x=0.");
    return TSIL_Infinity;
  }

  if (TSIL_CABS (1.0L - S/X) < 10.0L*TSIL_TOL) {
    TSIL_Warn("V0x0x", "V(0,x,0,x) is undefined for s=x.");
    return TSIL_Infinity;
  }

  lnbarX = TSIL_LOG(X/QQ);

  if (TSIL_CABS(S) < TSIL_TOL)
    return (2.0L - lnbarX - Zeta2)/X;

  S = AddIeps(S);

  log1mSoX = TSIL_CLOG(1.0L - S/X);

  return (-(Dilog(S/X) + Zeta2)/X  
	 + log1mSoX*((lnbarX -2.L+ 0.5L*log1mSoX)/S - 0.5L*log1mSoX/X));
}

/* ******************************************************************* */
/* special case of hep-ph/0307101 eq. (3.22)                           */

TSIL_COMPLEX V0x0y (TSIL_REAL X, TSIL_REAL Y, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_REAL lnbarY;

  if (TSIL_FABS(X) < TSIL_TOL) {
    TSIL_Warn("V0x0y", "V(0,x,0,y) is undefined for x=0.");
    return TSIL_Infinity;
  }

  if (TSIL_FABS(Y) < TSIL_TOL)
    return V0x00(X, S, QQ);

  if (TSIL_CABS(1-S/X) < 10.0*TSIL_TOL) {
    TSIL_Warn("V0x0y", "V(0,x,0,y) is undefined for s=x.");
    return TSIL_Infinity;
  }

  if (TSIL_FABS(X-Y) < TSIL_TOL)
    return V0x0x (X, S, QQ);

  lnbarY = TSIL_LOG(Y/QQ);

  return ((S*Y-X*X)*U0x0y(X,Y,S,QQ) + Y*(S-Y)*Tx00(Y,S,QQ) +
	  X*X -1.25L*S*X + 0.25L*S*Y + 2.L*X*Y - 2.L*Y*Y)/(X*(X - S)*(X - Y))
         + (2.0L*S00x(Y,S,QQ) - I20xy(X,Y,QQ) - Y*lnbarY)/(X*(X - S))
         + (X - Y*lnbarY + Y)*B0x(X,S,QQ)/(X*(X - Y));
}

/* ******************************************************************* */
/* special case of hep-ph/0307101 eq. (3.22)                           */

TSIL_COMPLEX V0xyz (TSIL_REAL X, TSIL_REAL Y, TSIL_REAL Z, TSIL_COMPLEX S, TSIL_REAL QQ)
{
  TSIL_REAL    temp, AX, AY, AZ, X2, X3, Y2, Y3, Z2, Z3;
  TSIL_COMPLEX denom;

  if (TSIL_FABS(Y) < TSIL_FABS(Z)) {temp=Z; Z=Y; Y=temp;}
 
  if (TSIL_FABS(X) < TSIL_TOL) {
    TSIL_Warn("V0xyz", "V(0,x,y,z) is undefined for x=0.");
    return TSIL_Infinity;
  }

  if (TSIL_FABS(Y) < TSIL_TOL)
    return V0x00 (X, S, QQ);
  
  if (TSIL_FABS(Z) < TSIL_TOL)
    return V0x0y (X, Y, S, QQ);

  if (TSIL_CABS(1.0L - S/X) < TSIL_TOL) {
    TSIL_Warn("V0xyz", "V(0,x,y,z) is undefined for s=x.");
    return TSIL_Infinity;
  }

  AX = A(X,QQ);
  AY = A(Y,QQ);
  AZ = A(Z,QQ);
  denom = X*(X-S)*Delta(X,Y,Z);
  X2 = X*X;
  X3 = X*X2;
  Y2 = Y*Y;
  Y3 = Y*Y2;
  Z2 = Z*Z;
  Z3 = Z*Z2;

  /* The special case below was added 11/14/04, but is not needed
     because it will have already been trapped setV. */
  /*
  if (TSIL_FABS(1.0L - (Y + Z + 2.0L*TSIL_SQRT(Y*Z))/X) < TSIL_TOL)
    return (4.L*S0xy(Y,Z,S,QQ) 
          + Tx0y(Z,Y,S,QQ)*((4.5L*Z + 1.5L*Y) + ((Z - Y)*(0.5L*S - Z + Y))/X)
          + Tx0y(Y,Z,S,QQ)*((4.5L*Y + 1.5L*Z) + ((Y - Z)*(0.5L*S - Y + Z))/X)
          + B0x(X,S,QQ)*((X - 2.0L*S) + (X*X + S*(Y-Z))*AY/(2.0L*X*Y) +
                              (X*X + S*(Z-Y))*AZ/(2.0L*X*Z))
          -AX + (0.5L*AX*AY*(X + Y - Z))/(X*Y) 
          + (0.5L*AX*AZ*(X - Y +Z))/(X*Z) 
          + (0.5L*AZ*(X*Y - Y*Y - 5.L*X*Z + 2.L*Y*Z - Z*Z))/(X*Z) 
          + (0.5L*AY*(-5.L*X*Y - Y*Y + X*Z + 2.L*Y*Z - Z*Z))/(X*Y) 
          + (0.5L*(-5.L*S*X + X*X + 10.L*X*Y - Y*Y + 10.L*X*Z + 
             2.L*Y*Z - Z*Z))/X 
          + (0.5L*AY*AZ*(-X*X + Y*Y - 2.L*Y*Z + Z*Z))/(X*Y*Z))/(X*(X-S));
  */

  return ((-AY*X2 - AZ*X2 - 1.25*S*X2 + X3 + 2.*AY*X*Y + 2.*AZ*X*Y +
	   1.5*S*X*Y - AY*Y2 - AZ*Y2 - 0.25*S*Y2 - 2.*X*Y2 + Y3 + 2.*AY*X*Z +
	   2.*AZ*X*Z + 1.5*S*X*Z + 2.*AY*Y*Z + 2.*AZ*Y*Z + 0.5*S*Y*Z 
           - 4.*X*Y*Z -
	   Y2*Z - AY*Z2 - AZ*Z2 - 0.25*S*Z2 - 2.*X*Z2 - Y*Z2 + Z3 +
	   (AY*S*X + AZ*S*X - AY*X2 - AZ*X2 - S*X2 + X3 - AY*S*Y +
	    AZ*S*Y + AY*X*Y - AZ*X*Y + S*X*Y - X2*Y + AY*S*Z - AZ*S*Z -
	    AY*X*Z + AZ*X*Z + S*X*Z - X2*Z)*B0x(X,S,QQ) +
	   (-X2 + 2.*X*Y - Y2 + 2.*X*Z + 2.*Y*Z - Z2)*I2(X,Y,Z,QQ) +
	   (2.*X2 - 4.*X*Y + 2.*Y2 - 4.*X*Z - 4.*Y*Z + 2.*Z2)*S0xy(Y,Z,S,QQ) +
	   (S*X*Y - S*Y2 - X*Y2 + Y3 + S*Y*Z - 3.*X*Y*Z - 2.*Y2*Z + Y*Z2)*
	   Tx0y(Y,Z,S,QQ) + (S*X*Z + S*Y*Z - 3.*X*Y*Z + Y2*Z - S*Z2 - X*Z2 -
			     2.*Y*Z2 + Z3)*Tx0y(Z,Y,S,QQ) +
	   (-X3 + S*X*Y + X2*Y - S*Y2 + S*X*Z + X2*Z + 2.L*S*Y*Z - S*Z2)*
	   U0xyz(X,Y,Z,S,QQ) )/denom);
}
