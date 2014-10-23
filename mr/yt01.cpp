#include <tt.hpp>
std::complex<long double> tt::my01()
{     
      
      
    std::complex<long double> myt[4], mytret;

    myt[1]=Tsil::A(MMt,mu2);
    myt[2]=pow(MMt,-1);
   myt[3]=myt[1]*myt[2];
   myt[3]= - 1./3. + myt[3];

      mytret = 4*myt[3];
      return mytret;
}
