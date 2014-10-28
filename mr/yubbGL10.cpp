#include <bb.hpp>
std::complex<long double>
bb::mygl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[9], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMW,-1);
    aryubbGL[4]=Tsil::A(MMH,mu2);
    aryubbGL[5]=Tsil::A(MMt,mu2);
    aryubbGL[6]=pow(MMH,-1);
   aryubbGL[7]=3*aryubbGL[4] + 1./2.*MMt;
   aryubbGL[8]= - MMt*aryubbGL[6];
   aryubbGL[8]=1./8. + aryubbGL[8];
   aryubbGL[8]=aryubbGL[5]*aryubbGL[8];
   aryubbGL[7]=1./8.*aryubbGL[7] + 3*aryubbGL[8];

      yubbGLret = aryubbGL[7]*aryubbGL[3]*pow(aryubbGL[2],2)*
      aryubbGL[1];
      return yubbGLret;
}
