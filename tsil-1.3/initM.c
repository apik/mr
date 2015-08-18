/* Routines for setup and initial values of M-type functions */

#include "internal.h"

/* ***************************************************************** */

void TSIL_Case0 (TSIL_MTYPE *foo, TSIL_REAL qq)
{
  TSIL_REAL x, y, z, u, v;
  TSIL_REAL Ax, Ay, Az, Au, Av;
  TSIL_REAL Aosqx, Aosqy, Aosqz, Aosqu;
  TSIL_REAL tempcoM, sqrtx, sqrty, sqrtz, sqrtu;

  foo->extramassdim1 = 0;
  foo->extramassdim2 = 0;

  /* For convenience */
  x = foo->arg[0];
  y = foo->arg[1];
  z = foo->arg[2];
  u = foo->arg[3];
  v = foo->arg[4];

  Ax = TSIL_A(x, qq);
  Ay = TSIL_A(y, qq);
  Az = TSIL_A(z, qq);
  Au = TSIL_A(u, qq);
  Av = TSIL_A(v, qq);

  sqrtx = TSIL_SQRT(x);
  sqrty = TSIL_SQRT(y);
  sqrtz = TSIL_SQRT(z);
  sqrtu = TSIL_SQRT(u);

  Aosqx = TSIL_Alpha(x, qq);
  Aosqy = TSIL_Alpha(y, qq);
  Aosqz = TSIL_Alpha(z, qq);
  Aosqu = TSIL_Alpha(u, qq);

  foo->adenom[2] = v;
  foo->adenom[1] = -u*v + v*v - u*x - v*x + u*y - v*y - v*z + x*z - y*z;
  foo->adenom[0] = u*u*x - u*v*x + u*x*x - u*x*y + v*x*y + u*v*z - u*x*z - 
            u*y*z - v*y*z - x*y*z + y*y*z + y*z*z;

  foo->aMU[0]   = 0.L;
  foo->aMUs[0]  = 1.L;
  foo->bMU[0][0] = (z*(y - x - v) + sqrtz*sqrtx*(u - z - v))/2;
  foo->bMU[0][1] = (z*(y - x - v) - sqrtz*sqrtx*(u - z - v))/2;

  foo->aMU[1]   = 0.L;
  foo->aMUs[1]  = 1.L;
  foo->bMU[1][0] = (u*(x - y - v) + sqrtu*sqrty*(z - u - v))/2;
  foo->bMU[1][1] = (u*(x - y - v) - sqrtu*sqrty*(z - u - v))/2;

  foo->aMU[2]   = 0.L; 
  foo->aMUs[2]  = 1.L;
  foo->bMU[2][0] = (x*(u - z - v) + sqrtz*sqrtx*(y - x - v))/2;
  foo->bMU[2][1] = (x*(u - z - v) - sqrtz*sqrtx*(y - x - v))/2;

  foo->aMU[3]   = 0.L; 
  foo->aMUs[3]  = 1.L;
  foo->bMU[3][0] = (y*(z - u - v) + sqrtu*sqrty*(x - y - v))/2;
  foo->bMU[3][1] = (y*(z - u - v) - sqrtu*sqrty*(x - y - v))/2;

  foo->aMS   = 2*v;
  foo->bMS[0] = -2*(foo->bMU[0][0] + foo->bMU[2][0]);
  foo->bMS[1] = -2*(foo->bMU[0][1] + foo->bMU[2][1]);
  foo->bMS[2] = -2*(foo->bMU[1][0] + foo->bMU[3][0]);
  foo->bMS[3] = -2*(foo->bMU[1][1] + foo->bMU[3][1]);
  foo->cMSconst = -0.5L*(TSIL_I2(x,y,v,qq) + TSIL_I2(z,u,v,qq));

  foo->aMT[0]   = x*(2*v + z - u);
  foo->bMT[0][0] = -3*x*foo->bMU[0][0] -(2*x + z)*foo->bMU[2][0];
  foo->bMT[0][1] = -3*x*foo->bMU[0][1] -(2*x + z)*foo->bMU[2][1];
  foo->bMT[0][2] = x*foo->bMS[2]/2;
  foo->bMT[0][3] = x*foo->bMS[3]/2;

  foo->aMT[1]   = y*(2*v + u - z);
  foo->bMT[1][0] = y*foo->bMS[0]/2;
  foo->bMT[1][1] = y*foo->bMS[1]/2;
  foo->bMT[1][2] = -3*y*foo->bMU[1][0] -(2*y + u)*foo->bMU[3][0];
  foo->bMT[1][3] = -3*y*foo->bMU[1][1] -(2*y + u)*foo->bMU[3][1];

  foo->aMT[2]   = z*(2*v + x - y);
  foo->bMT[2][0] = -3*z*foo->bMU[2][0] -(2*z + x)*foo->bMU[0][0];
  foo->bMT[2][1] = -3*z*foo->bMU[2][1] -(2*z + x)*foo->bMU[0][1];
  foo->bMT[2][2] = z*foo->bMS[2]/2;
  foo->bMT[2][3] = z*foo->bMS[3]/2;

  foo->aMT[3]   = u*(2*v + y - x);
  foo->bMT[3][0] = u*foo->bMS[0]/2;
  foo->bMT[3][1] = u*foo->bMS[1]/2;
  foo->bMT[3][2] = -3*u*foo->bMU[3][0] -(2*u + y)*foo->bMU[1][0];
  foo->bMT[3][3] = -3*u*foo->bMU[3][1] -(2*u + y)*foo->bMU[1][1];

  foo->aMT5s  = v;
  foo->aMT[4]   = v*v;
  foo->bMT[4][0] = v*foo->bMS[0]/2;
  foo->bMT[4][1] = v*foo->bMS[1]/2;
  foo->bMT[4][2] = v*foo->bMS[2]/2;
  foo->bMT[4][3] = v*foo->bMS[3]/2;

  foo->dMB  = 0.L;
  foo->dMBs = 1.L;

  foo->aMB[0]   = Av - v;
  foo->bMB[0][0] = (sqrty - Aosqy)*(sqrty*(z - u - v) + sqrtu*(x - y - v))/2 +
                (sqrtu - Aosqu)*(sqrtu*(x - y - v) + sqrty*(z - u - v))/2;
  foo->bMB[0][1] = (sqrty - Aosqy)*(sqrty*(z - u - v) - sqrtu*(x - y - v))/2 + 
                (sqrtu - Aosqu)*(sqrtu*(x - y - v) - sqrty*(z - u - v))/2;

  foo->aMB[1]   = Av - v;
  foo->bMB[1][0] = (sqrtz - Aosqz)*(sqrtz*(y - x - v) + sqrtx*(u - z - v))/2 + 
                (sqrtx - Aosqx)*(sqrtx*(u - z - v) + sqrtz*(y - x - v))/2;
  foo->bMB[1][1] = (sqrtz - Aosqz)*(sqrtz*(y - x - v) - sqrtx*(u - z - v))/2 + 
                (sqrtx - Aosqx)*(sqrtx*(u - z - v) - sqrtz*(y - x - v))/2;

  foo->aMs  = -3*v/2;
  tempcoM   = 2*v + x + y + z + u - Ax - Ay - Az - Au - 2*Av;
  foo->aM   = (x - y)*(u - z) -(x + y + z + u)*v/2 + foo->aMS*tempcoM/2;
  foo->bM[0] = (tempcoM*foo->bMS[0] + (3*x+z)*foo->bMU[0][0] + (3*z+x)*foo->bMU[2][0])/2;
  foo->bM[1] = (tempcoM*foo->bMS[1] + (3*x+z)*foo->bMU[0][1] + (3*z+x)*foo->bMU[2][1])/2;
  foo->bM[2] = (tempcoM*foo->bMS[2] + (3*y+u)*foo->bMU[1][0] + (3*u+y)*foo->bMU[3][0])/2;
  foo->bM[3] = (tempcoM*foo->bMS[3] + (3*y+u)*foo->bMU[1][1] + (3*u+y)*foo->bMU[3][1])/2;

  return;
}

/* ***************************************************************** */

void TSIL_Case1a (TSIL_MTYPE *foo, TSIL_REAL qq)
{
  TSIL_REAL x, y, z, u, v;
  TSIL_REAL Ax, Ay, Az, Au, Av;
  TSIL_REAL Aosqx, Aosqy, Aosqz, Aosqu, Aosqv;
  TSIL_REAL sqrtx, sqrty, sqrtz, sqrtu, sqrtv;

  foo->extramassdim1 = -3;
  foo->extramassdim2 = -1;

  /* For convenience */
  x = foo->arg[0];
  y = foo->arg[1];
  z = foo->arg[2];
  u = foo->arg[3];
  v = foo->arg[4];

  Ax = TSIL_A(x, qq);
  Ay = TSIL_A(y, qq);
  Az = TSIL_A(z, qq);
  Au = TSIL_A(u, qq);
  Av = TSIL_A(v, qq);

  Aosqx = TSIL_Alpha(x, qq);
  Aosqy = TSIL_Alpha(y, qq);
  Aosqz = TSIL_Alpha(z, qq);
  Aosqu = TSIL_Alpha(u, qq);
  Aosqv = TSIL_Alpha(v, qq);

  sqrtx = TSIL_SQRT(x);
  sqrty = TSIL_SQRT(y);
  sqrtz = TSIL_SQRT(z);
  sqrtu = TSIL_SQRT(u);
  sqrtv = TSIL_SQRT(v);

  foo->adenom[2] = 0.L;
  foo->adenom[1] = sqrtv;
  foo->adenom[0] = sqrtx*(y - u) + sqrty*(x - z);

  foo->aMU[0]   = 1.L;
  foo->aMUs[0]  = 0.L;
  foo->bMU[0][0] = -0.5L * sqrtz * foo->THxz;
  foo->bMU[0][1] =  0.5L * sqrtz * foo->PSxz;

  foo->aMU[1]   = 1.L;
  foo->aMUs[1]  = 0.L;
  foo->bMU[1][0] = -0.5L * sqrtu * foo->THyu;
  foo->bMU[1][1] =  0.5L * sqrtu * foo->PSyu;

  foo->aMU[2]   = 1.L;
  foo->aMUs[2]  = 0.L;
  foo->bMU[2][0] = -0.5L * sqrtx * foo->THxz;
  foo->bMU[2][1] = -0.5L * sqrtx * foo->PSxz;

  foo->aMU[3]   = 1.L; 
  foo->aMUs[3]  = 0.L;
  foo->bMU[3][0] = -0.5L * sqrty * foo->THyu;
  foo->bMU[3][1] = -0.5L * sqrty * foo->PSyu;

  foo->aMS   = 0.L;
  foo->bMS[0] = sqrtx + sqrtz;
  foo->bMS[1] = sqrtx - sqrtz;
  foo->bMS[2] = sqrty + sqrtu;
  foo->bMS[3] = sqrty - sqrtu;
  foo->cMSconst = -0.5L*(TSIL_I2(x,y,v,qq) + TSIL_I2(z,u,v,qq));

  foo->aMT[0]   = 0.L;
  foo->bMT[0][0] = sqrtx*(x + z/2 + 3*sqrtx*sqrtz/2);
  foo->bMT[0][1] = sqrtx*(x + z/2 - 3*sqrtx*sqrtz/2);
  foo->bMT[0][2] = x*(sqrty + sqrtu)/2;
  foo->bMT[0][3] = x*(sqrty - sqrtu)/2;

  foo->aMT[1]   = 0.L;
  foo->bMT[1][0] = y*(sqrtx+sqrtz)/2;
  foo->bMT[1][1] = y*(sqrtx-sqrtz)/2;
  foo->bMT[1][2] = sqrty*(y + u/2 + 3*sqrty*sqrtu/2);
  foo->bMT[1][3] = sqrty*(y + u/2 - 3*sqrty*sqrtu/2);

  foo->aMT[2]   = 0.L;
  foo->bMT[2][0] = sqrtz*(z + x/2 + 3*sqrtz*sqrtx/2);
  foo->bMT[2][1] = sqrtz*(-z - x/2 + 3*sqrtz*sqrtx/2);
  foo->bMT[2][2] = z*(sqrty + sqrtu)/2;
  foo->bMT[2][3] = z*(sqrty - sqrtu)/2;

  foo->aMT[3]   = 0.L;
  foo->bMT[3][0] = u*(sqrtx+sqrtz)/2;
  foo->bMT[3][1] = u*(sqrtx-sqrtz)/2;
  foo->bMT[3][2] = sqrtu*(u + y/2 + 3*sqrty*sqrtu/2);
  foo->bMT[3][3] = sqrtu*(-u - y/2 + 3*sqrty*sqrtu/2);

  foo->aMT5s  = 0.L;
  foo->aMT[4]   = sqrtv;
  foo->bMT[4][0] = sqrtv*(sqrtx + sqrtz)*(sqrty + sqrtx)/2;
  foo->bMT[4][1] = sqrtv*(sqrtx - sqrtz)*(sqrty + sqrtx)/2;
  foo->bMT[4][2] = sqrtv*(sqrty + sqrtu)*(sqrty + sqrtx)/2;
  foo->bMT[4][3] = sqrtv*(sqrty - sqrtu)*(sqrty + sqrtx)/2;

  foo->dMB  = 1.L;
  foo->dMBs = 0.L;

  foo->aMB[0]   = Aosqv - sqrtv;
  foo->bMB[0][0] = foo->THyu*(Aosqy + Aosqu - sqrty - sqrtu)/2;
  foo->bMB[0][1] = foo->PSyu*(Aosqy - Aosqu - sqrty + sqrtu)/2;

  foo->aMB[1]   = Aosqv - sqrtv;
  foo->bMB[1][0] = foo->THxz*(Aosqx + Aosqz - sqrtx - sqrtz)/2;
  foo->bMB[1][1] = foo->PSxz*(Aosqx - Aosqz - sqrtx + sqrtz)/2;

  foo->aM   = -3.L*sqrtv/2.L;
  foo->aMs  = 0.L;
  foo->bM[0] = (sqrtx + sqrtz)*(3.L*y + z/2.L + u + 5.L*x/2.L + (4.L*sqrty - sqrtz)*sqrtx -
			       Ay - 2.L*Av - Az - Au - Ax)/2.L;
  foo->bM[1] = (sqrtx - sqrtz)*(3.L*y + z/2.L + u + 5.L*x/2.L + (4.L*sqrty + sqrtz)*sqrtx -
			       Ay - 2.L*Av - Az - Au - Ax)/2.L;
  foo->bM[2] = (sqrty + sqrtu)*(5.L*y/2.L + z + u/2.L + 3.L*x + (4.L*sqrtx - sqrtu)*sqrty -
			       Ay - 2.L*Av - Az - Au - Ax)/2.L;
  foo->bM[3] = (sqrty - sqrtu)*(5.L*y/2.L + z + u/2.L + 3.L*x + (4.L*sqrtx + sqrtu)*sqrty -
			       Ay - 2.L*Av - Az - Au - Ax)/2.L;

  return;
}

/* ***************************************************************** */

void TSIL_Case1b (TSIL_MTYPE *foo, TSIL_REAL qq)
{
  TSIL_REAL x, y, z, u, v;
  TSIL_REAL Ax, Ay, Az, Au, Av;
  TSIL_REAL Aosqx, Aosqy, Aosqz, Aosqu, Aosqv;
  TSIL_REAL sqrtx, sqrty, sqrtz, sqrtu, sqrtv;

  foo->extramassdim1 = -3;
  foo->extramassdim2 = -1;

  /* For convenience */
  x = foo->arg[0];
  y = foo->arg[1];
  z = foo->arg[2];
  u = foo->arg[3];
  v = foo->arg[4];

  Ax = TSIL_A(x, qq);
  Ay = TSIL_A(y, qq);
  Az = TSIL_A(z, qq);
  Au = TSIL_A(u, qq);
  Av = TSIL_A(v, qq);

  Aosqx = TSIL_Alpha(x, qq);
  Aosqy = TSIL_Alpha(y, qq);
  Aosqz = TSIL_Alpha(z, qq);
  Aosqu = TSIL_Alpha(u, qq);
  Aosqv = TSIL_Alpha(v, qq);

  sqrtx = TSIL_SQRT(x);
  sqrty = TSIL_SQRT(y);
  sqrtz = TSIL_SQRT(z);
  sqrtu = TSIL_SQRT(u);
  sqrtv = TSIL_SQRT(v);

  foo->adenom[2] = 0.L;
  foo->adenom[1] = sqrtv;
  foo->adenom[0] = sqrtz*(u - y) + sqrtu*(z-x);

  foo->aMU[0] = 1.L;
  foo->aMUs[0] = 0.L;
  foo->bMU[0][0] = -sqrtz * foo->THxz/2;
  foo->bMU[0][1] = -sqrtz * foo->PSxz/2;

  foo->aMU[1] = 1.L;
  foo->aMUs[1] = 0.L;
  foo->bMU[1][0] = -sqrtu * foo->THyu/2;
  foo->bMU[1][1] = -sqrtu * foo->PSyu/2;

  foo->aMU[2] = 1.L;
  foo->aMUs[2] = 0.L;
  foo->bMU[2][0] = -sqrtx * foo->THxz/2;
  foo->bMU[2][1] = sqrtx * foo->PSxz/2;

  foo->aMU[3] = 1.L;
  foo->aMUs[3] = 0.L;
  foo->bMU[3][0] = -sqrty * foo->THyu/2;
  foo->bMU[3][1] = sqrty * foo->PSyu/2;

  foo->aMS = 0.L;
  foo->bMS[0] = sqrtz + sqrtx;
  foo->bMS[1] = sqrtz - sqrtx;
  foo->bMS[2] = sqrtu + sqrty;
  foo->bMS[3] = sqrtu - sqrty;
  foo->cMSconst = -0.5L*(TSIL_I2(x,y,v,qq) + TSIL_I2(z,u,v,qq));

  foo->aMT[0] = 0.L;
  foo->bMT[0][0] = sqrtx*(x + z/2 + 3*sqrtx*sqrtz/2);
  foo->bMT[0][1] = sqrtx*(-x - z/2 + 3*sqrtx*sqrtz/2);
  foo->bMT[0][2] = x*(sqrtu + sqrty)/2;
  foo->bMT[0][3] = x*(sqrtu - sqrty)/2;

  foo->aMT[1] = 0.L;
  foo->bMT[1][0] = y*(sqrtz + sqrtx)/2;
  foo->bMT[1][1] = y*(sqrtz - sqrtx)/2;
  foo->bMT[1][2] = sqrty*(y + u/2 + 3*sqrtu*sqrty/2);
  foo->bMT[1][3] = sqrty*(-y - u/2 + 3*sqrtu*sqrty/2);

  foo->aMT[2] = 0.L;
  foo->bMT[2][0] = sqrtz*(z + x/2 + 3*sqrtz*sqrtx/2);
  foo->bMT[2][1] = sqrtz*(z + x/2 - 3*sqrtz*sqrtx/2);
  foo->bMT[2][2] = z*(sqrtu + sqrty)/2;
  foo->bMT[2][3] = z*(sqrtu - sqrty)/2;

  foo->aMT[3] = 0.L;
  foo->bMT[3][0] = u*(sqrtz+sqrtx)/2;
  foo->bMT[3][1] = u*(sqrtz-sqrtx)/2;
  foo->bMT[3][2] = sqrtu*(u + y/2 + 3*sqrtu*sqrty/2);
  foo->bMT[3][3] = sqrtu*(u + y/2 - 3*sqrtu*sqrty/2);

  foo->aMT5s = 0.L;
  foo->aMT[4] = sqrtv;
  foo->bMT[4][0] = sqrtv*(sqrtz + sqrtx)*(sqrtu + sqrtz)/2;
  foo->bMT[4][1] = sqrtv*(sqrtz - sqrtx)*(sqrtu + sqrtz)/2;
  foo->bMT[4][2] = sqrtv*(sqrtu + sqrty)*(sqrtu + sqrtz)/2;
  foo->bMT[4][3] = sqrtv*(sqrtu - sqrty)*(sqrtu + sqrtz)/2;

  foo->dMB = 1.L;
  foo->dMBs = 0.L;

  foo->aMB[0] = Aosqv - sqrtv;
  foo->bMB[0][0] = foo->THyu*(Aosqu + Aosqy - sqrtu - sqrty)/2;
  foo->bMB[0][1] = foo->PSyu*(Aosqu - Aosqy - sqrtu + sqrty)/2;

  foo->aMB[1] = Aosqv - sqrtv;
  foo->bMB[1][0] = foo->THxz*(Aosqz + Aosqx - sqrtz - sqrtx)/2;
  foo->bMB[1][1] = foo->PSxz*(Aosqz - Aosqx - sqrtz + sqrtx)/2;

  foo->aM = -3*sqrtv/2;
  foo->aMs = 0.L;
  foo->bM[0] = (sqrtz + sqrtx)*(3*u + x/2 + y + 5*z/2 + (4*sqrtu - sqrtx)*sqrtz - 
			       Au - 2*Av - Ax - Ay - Az)/2;
  foo->bM[1] = (sqrtz - sqrtx)*(3*u + x/2 + y + 5*z/2 + (4*sqrtu + sqrtx)*sqrtz - 
			       Au - 2*Av - Ax - Ay - Az)/2;
  foo->bM[2] = (sqrtu + sqrty)*(5*u/2 + x + y/2 + 3*z + (4*sqrtz - sqrty)*sqrtu - 
			       Au - 2*Av - Ax - Ay - Az)/2;
  foo->bM[3] = (sqrtu - sqrty)*(5*u/2 + x + y/2 + 3*z + (4*sqrtz + sqrty)*sqrtu - 
			       Au - 2*Av - Ax - Ay - Az)/2;

  return;
}

/* ***************************************************************** */

void TSIL_Case2a (TSIL_MTYPE *foo, TSIL_REAL qq)
{
  TSIL_REAL x, y, z, u, v;
  TSIL_REAL Ay, Au, Av;
  TSIL_REAL Aosqy, Aosqu;
  TSIL_REAL sqrty, sqrtu;

  foo->extramassdim1 = -2;
  foo->extramassdim2 = 0;

  /* For convenience */
  x = foo->arg[0];
  y = foo->arg[1];
  z = foo->arg[2];
  u = foo->arg[3];
  v = foo->arg[4];

  Ay = TSIL_A(y, qq);
  Au = TSIL_A(u, qq);
  Av = TSIL_A(v, qq);

  Aosqy = TSIL_Alpha(y, qq);
  Aosqu = TSIL_Alpha(u, qq);

  sqrty = TSIL_SQRT(y);
  sqrtu = TSIL_SQRT(u);

  foo->adenom[2] = 0.L;
  foo->adenom[1] = v;
  foo->adenom[0] = (u-v)*(y-v);

  foo->aMU[0] = 0.L;
  foo->aMUs[0] = 0.L;
  foo->bMU[0][0] = 0.L;
  foo->bMU[0][1] = 0.L;

  foo->aMU[1] = 1.L;
  foo->aMUs[1] = 0.L;
  foo->bMU[1][0] = -sqrtu*(sqrtu + sqrty)*(v + sqrtu*sqrty)/2;
  foo->bMU[1][1] = sqrtu*(sqrtu - sqrty)*(-v + sqrtu*sqrty)/2;

  foo->aMU[2] = 0.L; 
  foo->aMUs[2] = 0.L;
  foo->bMU[2][0] = 0.L;
  foo->bMU[2][1] = 0.L;

  foo->aMU[3] = 1.L; 
  foo->aMUs[3] = 0.L;
  foo->bMU[3][0] = -sqrty*(sqrtu + sqrty)*(v + sqrtu*sqrty)/2;
  foo->bMU[3][1] = sqrty*(sqrtu - sqrty)*(v - sqrtu*sqrty)/2;

  foo->aMS = 0.L;
  foo->bMS[0] = 0.L;
  foo->bMS[1] = 0.L;
  foo->bMS[2] = v + sqrtu*sqrty;
  foo->bMS[3] = v - sqrtu*sqrty;
  foo->cMSconst = -0.5L*(TSIL_I2(x,y,v,qq) + TSIL_I2(z,u,v,qq));

  foo->aMT[0] = 0.L;
  foo->bMT[0][0] = 0.L;
  foo->bMT[0][1] = 0.L;
  foo->bMT[0][2] = 0.L;
  foo->bMT[0][3] = 0.L;

  foo->aMT[1] = 0.L;
  foo->bMT[1][0] = 0.L;
  foo->bMT[1][1] = 0.L;
  foo->bMT[1][2] = sqrty*(2*sqrty + sqrtu)*(v + sqrty*sqrtu)/2;
  foo->bMT[1][3] = sqrty*(2*sqrty - sqrtu)*(v - sqrty*sqrtu)/2;

  foo->aMT[2] = 0.L;
  foo->bMT[2][0] = 0.L;
  foo->bMT[2][1] = 0.L;
  foo->bMT[2][2] = 0.L;
  foo->bMT[2][3] = 0.L;

  foo->aMT[3] = 0.L;
  foo->bMT[3][0] = 0.L;
  foo->bMT[3][1] = 0.L;
  foo->bMT[3][2] = sqrtu*(2*sqrtu + sqrty)*(v + sqrty*sqrtu)/2;
  foo->bMT[3][3] = sqrtu*(2*sqrtu - sqrty)*(v - sqrty*sqrtu)/2;

  foo->aMT5s = 0.L;
  foo->aMT[4] = v;
  foo->bMT[4][0] = 0.L;
  foo->bMT[4][1] = 0.L;
  foo->bMT[4][2] = v*(v + sqrtu*sqrty)/2;
  foo->bMT[4][3] = v*(v - sqrtu*sqrty)/2;

  foo->dMB = 1.L;
  foo->dMBs = 0.L;

  foo->aMB[0] = Av - v;
  foo->bMB[0][0] = (sqrty + sqrtu)*(v + sqrtu*sqrty)*(Aosqy + Aosqu - sqrtu - sqrty)/2;
  foo->bMB[0][1] = (sqrty - sqrtu)*(v - sqrtu*sqrty)*(Aosqy - Aosqu + sqrtu - sqrty)/2;

  foo->aMB[1] = Av - v;
  foo->bMB[1][0] = 0.L;
  foo->bMB[1][1] = 0.L;

  foo->aM = -3.L*v/2.L;
  foo->aMs = 0.L;
  foo->bM[0] = 0.L;
  foo->bM[1] = 0.L;
  foo->bM[2] = (v + sqrty*sqrtu)*(u + y + 4*v - 2*Au - 2*Ay - 4*Av - 2*sqrtu*sqrty)/4;
  foo->bM[3] = (v - sqrty*sqrtu)*(u + y + 4*v - 2*Au - 2*Ay - 4*Av + 2*sqrtu*sqrty)/4;

  return;
}

/* ***************************************************************** */

void TSIL_Case2b (TSIL_MTYPE *foo, TSIL_REAL qq)
{
  TSIL_REAL x, y, z, u, v;
  TSIL_REAL Ax, Az, Av;
  TSIL_REAL Aosqx, Aosqz;
  TSIL_REAL sqrtx, sqrtz;

  foo->extramassdim1 = -2;
  foo->extramassdim2 = 0;

  /* For convenience */
  x = foo->arg[0];
  y = foo->arg[1];
  z = foo->arg[2];
  u = foo->arg[3];
  v = foo->arg[4];

  Ax = TSIL_A(x, qq);
  Az = TSIL_A(z, qq);
  Av = TSIL_A(v, qq);

  Aosqx = TSIL_Alpha(x, qq);
  Aosqz = TSIL_Alpha(z, qq);

  sqrtx = TSIL_SQRT(x);
  sqrtz = TSIL_SQRT(z);

  foo->adenom[2] = 0.L;
  foo->adenom[1] = v;
  foo->adenom[0] = (z-v)*(x-v);

  foo->aMU[0] = 1.L;
  foo->aMUs[0] = 0.L;
  foo->bMU[0][0] = -sqrtz*(sqrtz + sqrtx)*(v + sqrtz*sqrtx)/2;
  foo->bMU[0][1] = sqrtz*(sqrtz - sqrtx)*(-v + sqrtz*sqrtx)/2;

  foo->aMU[1] = 0.L;
  foo->aMUs[1] = 0.L;
  foo->bMU[1][0] = 0.L;
  foo->bMU[1][1] = 0.L;

  foo->aMU[2] = 1.L;
  foo->aMUs[2] = 0.L;
  foo->bMU[2][0] = -sqrtx*(sqrtz + sqrtx)*(v + sqrtz*sqrtx)/2;
  foo->bMU[2][1] = sqrtx*(sqrtz - sqrtx)*(v - sqrtz*sqrtx)/2;

  foo->aMU[3] = 0.L;
  foo->aMUs[3] = 0.L;
  foo->bMU[3][0] = 0.L;
  foo->bMU[3][1] = 0.L;

  foo->aMS = 0.L;
  foo->bMS[0] = v + sqrtz*sqrtx;
  foo->bMS[1] = v - sqrtz*sqrtx;
  foo->bMS[2] = 0.L;
  foo->bMS[3] = 0.L;
  foo->cMSconst = -0.5L*(TSIL_I2(x,y,v,qq) + TSIL_I2(z,u,v,qq));

  foo->aMT[0] = 0.L;
  foo->bMT[0][0] = sqrtx*(2*sqrtx + sqrtz)*(v + sqrtx*sqrtz)/2;
  foo->bMT[0][1] = sqrtx*(2*sqrtx - sqrtz)*(v - sqrtx*sqrtz)/2;
  foo->bMT[0][2] = 0.L;
  foo->bMT[0][3] = 0.L;

  foo->aMT[1] = 0.L;
  foo->bMT[1][0] = 0.L;
  foo->bMT[1][1] = 0.L;
  foo->bMT[1][2] = 0.L;
  foo->bMT[1][3] = 0.L;

  foo->aMT[2] = 0.L;
  foo->bMT[2][0] = sqrtz*(2*sqrtz + sqrtx)*(v + sqrtx*sqrtz)/2;
  foo->bMT[2][1] = sqrtz*(2*sqrtz - sqrtx)*(v - sqrtx*sqrtz)/2;
  foo->bMT[2][2] = 0.L;
  foo->bMT[2][3] = 0.L;

  foo->aMT[3] = 0.L;
  foo->bMT[3][0] = 0.L;
  foo->bMT[3][1] = 0.L;
  foo->bMT[3][2] = 0.L;
  foo->bMT[3][3] = 0.L;

  foo->aMT5s = 0.L;
  foo->aMT[4] = v;
  foo->bMT[4][0] = v*(v + sqrtz*sqrtx)/2;
  foo->bMT[4][1] = v*(v - sqrtz*sqrtx)/2;
  foo->bMT[4][2] = 0.L;
  foo->bMT[4][3] = 0.L;

  foo->dMB = 1.L;
  foo->dMBs = 0.L;

  foo->aMB[0] = Av - v;
  foo->bMB[0][0] = 0.L;
  foo->bMB[0][1] = 0.L;

  foo->aMB[1] = Av - v;
  foo->bMB[1][0] = (sqrtx + sqrtz)*(v + sqrtz*sqrtx)*(Aosqx + Aosqz - sqrtz - sqrtx)/2;
  foo->bMB[1][1] = (sqrtx - sqrtz)*(v - sqrtz*sqrtx)*(Aosqx - Aosqz + sqrtz - sqrtx)/2;

  foo->aM = -3*v/2;
  foo->aMs = 0.L;
  foo->bM[2] = 0.L;
  foo->bM[3] = 0.L;
  foo->bM[0] = (v + sqrtx*sqrtz)*(z + x + 4*v - 2*Az - 2*Ax - 4*Av - 2*sqrtz*sqrtx)/4;
  foo->bM[1] = (v - sqrtx*sqrtz)*(z + x + 4*v - 2*Az - 2*Ax - 4*Av + 2*sqrtz*sqrtx)/4;

  return;
}

/* ***************************************************************** */

void TSIL_Case3 (TSIL_MTYPE *foo, TSIL_REAL qq)
{
  TSIL_REAL x, y, z, u, v;
  TSIL_REAL Ax, Ay, Av;

  foo->extramassdim1 = -2;
  foo->extramassdim2 = 0;

  /* For convenience */
  x = foo->arg[0];
  y = foo->arg[1];
  z = foo->arg[2];
  u = foo->arg[3];
  v = foo->arg[4];

  Ax = TSIL_A(x, qq);
  Ay = TSIL_A(y, qq);
  Av = TSIL_A(v, qq);

  foo->adenom[2] = 0.L;
  foo->adenom[1] = v;
  foo->adenom[0] = x*x + y*y + v*v - 2*x*y - 2*x*v - 2*y*v;

  foo->aMU[0] = 1.L;
  foo->aMUs[0] = 0.L;
  foo->bMU[0][0] = 2*x*(y - x - v);
  foo->bMU[0][1] = 0.L;

  foo->aMU[1] = 1.L;
  foo->aMUs[1] = 0.L;
  foo->bMU[1][0] = 2*y*(x - y - v);
  foo->bMU[1][1] = 0.L;

  foo->aMU[2] = 0.L;
  foo->aMUs[2] = 0.L;
  foo->bMU[2][0] = 0.L;
  foo->bMU[2][1] = 0.L;

  foo->aMU[3] = 0.L;
  foo->aMUs[3] = 0.L;
  foo->bMU[3][0] = 0.L;
  foo->bMU[3][1] = 0.L;

  foo->aMS = 0.L;
  foo->bMS[0] = (v + x - y);
  foo->bMS[1] = 0.L;
  foo->bMS[2] = (v - x + y);
  foo->bMS[3] = 0.L;
  foo->cMSconst = -0.5L*(TSIL_I2(x,y,v,qq) + TSIL_I2(z,u,v,qq));

  foo->aMT[0] = 0.L;
  foo->bMT[0][0] = 3*x*(v + x - y);
  foo->bMT[0][1] = 0.L;
  foo->bMT[0][2] = x*(v - x + y);
  foo->bMT[0][3] = 0.L;

  foo->aMT[1] = 0.L;
  foo->bMT[1][0] = y*(v - y + x);
  foo->bMT[1][1] = 0.L;
  foo->bMT[1][2] = 3*y*(v - x + y);
  foo->bMT[1][3] = 0.L;

  foo->aMT[2] = 0.L;
  foo->bMT[2][0] = 0.L;
  foo->bMT[2][1] = 0.L;
  foo->bMT[2][2] = 0.L;
  foo->bMT[2][3] = 0.L;

  foo->aMT[3] = 0.L;
  foo->bMT[3][0] = 0.L;
  foo->bMT[3][1] = 0.L;
  foo->bMT[3][2] = 0.L;
  foo->bMT[3][3] = 0.L;

  foo->aMT5s = 0.L;
  foo->aMT[4] = v;
  foo->bMT[4][0] = v*(v + x - y)/2;
  foo->bMT[4][1] = 0.L;
  foo->bMT[4][2] = v*(v - x + y)/2;
  foo->bMT[4][3] = 0.L;

  foo->dMB = 1.L;
  foo->dMBs = 0.L;

  foo->aMB[0] = Av - v;
  foo->bMB[0][0] = 2*(v - x + y)*(Ay - y);
  foo->bMB[0][1] = 0.L;

  foo->aMB[1] = Av - v;
  foo->bMB[1][0] = 2*(v + x - y)*(Ax - x);
  foo->bMB[1][1] = 0.L;

  foo->aM = -3*v/2;
  foo->aMs = 0.L;
  foo->bM[0] = (v + x - y)*(v + y - Ax - Ay - Av);
  foo->bM[1] = 0.L;
  foo->bM[2] = (v - x + y)*(v + x - Ax - Ay - Av);
  foo->bM[3] = 0.L;

  return;
}

/* ***************************************************************** */

int TSIL_ConstructM (TSIL_MTYPE *m,
		TSIL_REAL x,
		TSIL_REAL y, 
		TSIL_REAL z,
		TSIL_REAL u,
		TSIL_REAL v, 
		TSIL_REAL qq)
{
  m->arg[0] = x;
  m->arg[1] = y;
  m->arg[2] = z;
  m->arg[3] = u;
  m->arg[4] = v;

  m->THxz = TSIL_Th2(x, z);
  m->THyu = TSIL_Th2(y, u);
  m->PSxz = TSIL_Ps2(x, z);
  m->PSyu = TSIL_Ps2(y, u);

  if (TSIL_FABS(v - TSIL_Th2(x, y)) < TSIL_TOL) {
    TSIL_Case1a (m, qq);
    return 1;
  }
  else if (TSIL_FABS(v - TSIL_Th2(z, u)) < TSIL_TOL) {
    TSIL_Case1b (m, qq);
    return 1;
  }
  else if (x < TSIL_TOL && z < TSIL_TOL) {
    TSIL_Case2a (m, qq);
    return 2;
  }
  else if (y < TSIL_TOL && u < TSIL_TOL) {
    TSIL_Case2b (m, qq);
    return 2;
  }
  else if (TSIL_FABS(x - z) < TSIL_TOL && TSIL_FABS(y - u) < TSIL_TOL) {
    TSIL_Case3 (m, qq);
    return 3;
  }
  else {
    TSIL_Case0 (m, qq);
    return 0;
  }
}

