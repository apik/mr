#include <tt.hpp>
std::complex<long double>
tt::mygl01(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[5], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=Tsil::A(MMt,mu2);
    aryuttGL[3]=pow(MMt,-1);
   aryuttGL[4]=aryuttGL[2]*aryuttGL[3];
   aryuttGL[4]= - 1./3. + aryuttGL[4];

      yuttGLret = 4*aryuttGL[4]*aryuttGL[1];
      return yuttGLret;
}
