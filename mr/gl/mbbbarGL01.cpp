#include <bb.hpp>
long double bb::xgl01(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbarGL[5], mbbbarGLret;

    armbbbarGL[1]=double(boson);
    armbbbarGL[2]=Tsil::A(mmb,mu2);
    armbbbarGL[3]=pow(mmb,-1);
   armbbbarGL[4]=armbbbarGL[3]*armbbbarGL[2];
   armbbbarGL[4]= - 1./3. + armbbbarGL[4];

      mbbbarGLret = 4*armbbbarGL[4]*armbbbarGL[1];
      return mbbbarGLret.real();
}
