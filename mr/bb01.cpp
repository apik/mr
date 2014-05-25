#include <bb.hpp>
std::complex<long double> bb::m01(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbb[3];

    mbb[1]=Tsil::A(MMb,mu2);
    mbb[2]=pow(MMb,-1);
   mbb[3]=mbb[1]*mbb[2];
   mbb[3]= - 1./3. + mbb[3];

      return 4*mbb[3];
}
