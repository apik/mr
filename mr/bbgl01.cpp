#include <bb.hpp>
std::complex<long double> bb::mgl01(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbbgl[3];

    mbbgl[1]=log(pow(mu2,-1)*MMb);
   mbbgl[2]= - 4./3. + mbbgl[1];

      return 4*mbbgl[2];
}
