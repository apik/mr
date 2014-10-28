#include <tt.hpp>
#include "CRunDec.h"
std::complex<long double> tt::m03(size_t nL)
{     
  return std::complex<long double>(64.*(CRunDec::fMsFromOs3(sqrt(mu2),sqrt(MMt),nL) + CRunDec::fZmM(nL)), 0);
}
