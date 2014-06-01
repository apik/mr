#include <tt.hpp>
#include "CRunDec.h"
std::complex<long double> tt::m03(size_t nL)
{     
  long double ret = 64.*CRunDec::fMsFromOs3(sqrt(mu2),sqrt(MMt),nL);
  std::cout << "3-loop QCD: " << ret << std::endl;
  return std::complex<long double>(ret, 0);
}
