#include <tt.hpp>
std::complex<long double> tt::mygl01(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[3];

    mytgl[1]=log(pow(mu2,-1)*MMt);
   mytgl[2]= - 4./3. + mytgl[1];

      return 4*mytgl[2];
}
