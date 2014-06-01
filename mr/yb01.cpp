#include <bb.hpp>
std::complex<long double> bb::my01(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myb[4];

    myb[1]=Tsil::A(MMb,mu2);
    myb[2]=pow(MMb,-1);
   myb[3]=myb[1]*myb[2];
   myb[3]= - 1./3. + myb[3];

      return 4*myb[3];
}
