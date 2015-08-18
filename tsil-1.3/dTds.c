/* Derivative of T function with respect to s. */

#include "internal.h"

/* **************************************************************** */

TSIL_COMPLEX TSIL_dTds (TSIL_TTYPE T, TSIL_COMPLEX s)
{
  int i;
  TSIL_COMPLEX result = 0.0L;
  
  for (i=0; i<4; i++)
  {
    result += ( *(T.sval) * T.cTS_num[i] 
              + *(T.tval[0]) * T.cTT1_num[i] 
              + *(T.tval[1]) * T.cTT2_num[i] 
              + *(T.tval[2]) * T.cTT3_num[i]  
              + T.cT_num[i])/(s - T.T_den[i]);
  }

  result /= s;
  result += T.cTs_num/(s - T.T_den[0]);

  return result;
}
