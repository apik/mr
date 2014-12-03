#include <tt.hpp>
long double tt::xgl01(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[5], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=Tsil::A(MMt,mu2);
    armttbarGL[3]=pow(MMt,-1);
   armttbarGL[4]=armttbarGL[3]*armttbarGL[2];
   armttbarGL[4]= - 1./3. + armttbarGL[4];

      mttbarGLret = 4*armttbarGL[4]*armttbarGL[1];
      return mttbarGLret.real();
}
