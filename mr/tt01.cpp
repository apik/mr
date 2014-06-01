#include <tt.hpp>
std::complex<long double> tt::m01()
{     
      
      
    std::complex<long double> mtt[4];

    mtt[1]=Tsil::A(MMt,mu2);
    mtt[2]=pow(MMt,-1);
   mtt[3]=mtt[1]*mtt[2];
   mtt[3]= - 1./3. + mtt[3];

      return 4*mtt[3];
}
