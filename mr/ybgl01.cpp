#include <bb.hpp>
std::complex<long double> bb::mygl01(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mybgl[3];

    mybgl[1]=log(pow(mu2,-1)*MMb);
   mybgl[2]= - 4./3. + mybgl[1];

      return 4*mybgl[2];
}
