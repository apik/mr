#include <bb.hpp>
std::complex<long double>
bb::mygl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[6], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMW,-1);
    aryubbGL[4]=Tsil::A(MMt,mu2);
   aryubbGL[5]=1./2.*MMH - 3*aryubbGL[4] - 5./2.*MMt;

      yubbGLret = 1./8.*aryubbGL[5]*aryubbGL[3]*pow(aryubbGL[2],2)*
      aryubbGL[1];
      return yubbGLret;
}
