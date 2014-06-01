#include <tt.hpp>
#include "CRunDec.h"
std::complex<long double> tt::m02(size_t nL)
{     
  return std::complex<long double>(16.*CRunDec::fMsFromOs2(sqrt(mu2),sqrt(MMt),nL));
}
