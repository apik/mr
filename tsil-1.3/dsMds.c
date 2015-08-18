/* Derivative of s*M with respect to s. */

#include "internal.h"

/* ******************************************************************* */

TSIL_COMPLEX TSIL_dsMds (TSIL_MTYPE m, TSIL_COMPLEX s)
{
  TSIL_COMPLEX res;
  int i;
  TSIL_COMPLEX Mdenom;
  TSIL_COMPLEX cMU[4], cMT[5], cMB[2], cMS, cM;
  TSIL_COMPLEX invsminusTHxz, invsminusPSxz, invsminusTHyu, invsminusPSyu;

  invsminusTHxz = 1.L/(s - m.THxz);
  invsminusPSxz = 1.L/(s - m.PSxz);
  invsminusTHyu = 1.L/(s - m.THyu);
  invsminusPSyu = 1.L/(s - m.PSyu);

  cMU[0] = (m.aMU[0] + s*m.aMUs[0])*
           (m.bMU[0][0]*invsminusTHxz + m.bMU[0][1]*invsminusPSxz);
  cMU[1] = (m.aMU[1] + s*m.aMUs[1])*
           (m.bMU[1][0]*invsminusTHyu + m.bMU[1][1]*invsminusPSyu);
  cMU[2] = (m.aMU[2] + s*m.aMUs[2])*
           (m.bMU[2][0]*invsminusTHxz + m.bMU[2][1]*invsminusPSxz);
  cMU[3] = (m.aMU[3] + s*m.aMUs[3])*
           (m.bMU[3][0]*invsminusTHyu + m.bMU[3][1]*invsminusPSyu);

  for (i=0; i<5; i++)
    cMT[i] = m.aMT[i] + m.bMT[i][0]*invsminusTHxz + m.bMT[i][1]*invsminusPSxz 
                      + m.bMT[i][2]*invsminusTHyu + m.bMT[i][3]*invsminusPSyu;

  cMT[4] += s*m.aMT5s; 

  cMS = m.aMS + m.bMS[0]*invsminusTHxz + m.bMS[1]*invsminusPSxz
              + m.bMS[2]*invsminusTHyu + m.bMS[3]*invsminusPSyu;

  cMB[0] = (m.aMB[0] + m.bMB[0][0]*invsminusTHyu + m.bMB[0][1]*invsminusPSyu);

  cMB[1] = (m.aMB[1] + m.bMB[1][0]*invsminusTHxz + m.bMB[1][1]*invsminusPSxz);

  cM = m.aM + s*m.aMs + m.bM[0]*invsminusTHxz + m.bM[1]*invsminusPSxz 
                      + m.bM[2]*invsminusTHyu + m.bM[3]*invsminusPSyu;

  /* Assemble result: */
  res =  cMT[4] * (*(m.tval[4]) + *(m.tval[5]))
        + (m.dMB + s*m.dMBs)*(cMB[0] * *(m.bval[0]) + cMB[1] * *(m.bval[1]))
        + cM + cMS  * (*(m.sval[0]) + *(m.sval[1]) 
                      + 0.5L*s * *(m.bval[0]) * *(m.bval[1]) + m.cMSconst);

  for (i=0; i<4; i++) {
    res += cMU[i] * *(m.uval[i]);
    res += cMT[i] * *(m.tval[i]);
  }

  Mdenom = m.adenom[0] + m.adenom[1]*s + m.adenom[2]*s*s;

  return res/Mdenom;
}
