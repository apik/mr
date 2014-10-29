#include <bb.hpp>
std::complex<long double>
bb::mygl01(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[5], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=Tsil::A(mmb,mu2);
    aryubbGL[3]=pow(mmb,-1);
   aryubbGL[4]=aryubbGL[3]*aryubbGL[2];
   aryubbGL[4]= - 1./3. + aryubbGL[4];

      yubbGLret = 4*aryubbGL[4]*aryubbGL[1];
      return yubbGLret;
}
