#include <tt.hpp>
std::complex<long double> tt::mygl01(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[4];

    mytgl[1]=Tsil::A(MMt,mu2);
    mytgl[2]=pow(MMt,-1);
   mytgl[3]=mytgl[1]*mytgl[2];
   mytgl[3]= - 1./3. + mytgl[3];

      return 4*mytgl[3];
}
