#include <tt.hpp>
std::complex<long double> tt::mgl01(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mttgl[4];

    mttgl[1]=Tsil::A(MMt,mu2);
    mttgl[2]=pow(MMt,-1);
   mttgl[3]=mttgl[1]*mttgl[2];
   mttgl[3]= - 1./3. + mttgl[3];

      return 4*mttgl[3];
}
