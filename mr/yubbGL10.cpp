#include <bb.hpp>
long double bb<OS>::ygl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[6], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMW,-1);
    aryubbGL[4]=Tsil::A(MMt,mu2);
   aryubbGL[5]=MMH - 5*MMt;
   aryubbGL[5]=1./2.*aryubbGL[5] - 3*aryubbGL[4];

      yubbGLret = 1./8.*aryubbGL[5]*aryubbGL[3]*pow(aryubbGL[2],2)*
      aryubbGL[1];
      return yubbGLret.real();
}
