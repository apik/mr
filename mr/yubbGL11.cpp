#include <bb.hpp>
std::complex<long double>
bb::mygl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[17], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMH,-1);
    aryubbGL[4]=pow(MMW,-1);
    aryubbGL[5]=Tsil::I2(0,0,MMt,mu2);
    aryubbGL[6]=Tsil::A(MMH,mu2);
    aryubbGL[7]=Tsil::A(MMb,mu2);
    aryubbGL[8]=pow(MMb,-1);
    aryubbGL[9]=Tsil::A(MMt,mu2);
    aryubbGL[10]=pow(MMt,-1);
    aryubbGL[11]=Tsil::Aeps(MMt,mu2);
    aryubbGL[12]=Tsil::Aeps(MMb,mu2);
   aryubbGL[13]=MMt*aryubbGL[3];
   aryubbGL[14]= - MMt*aryubbGL[3];
   aryubbGL[14]=1./2. + 4*aryubbGL[14];
   aryubbGL[14]=aryubbGL[7]*aryubbGL[8]*aryubbGL[14];
   aryubbGL[15]= - 1./2.*aryubbGL[10] - 24*aryubbGL[3];
   aryubbGL[15]=aryubbGL[9]*aryubbGL[15];
   aryubbGL[14]=aryubbGL[15] + 3*aryubbGL[14] + 13./2. + 4*aryubbGL[13]
   ;
   aryubbGL[14]=aryubbGL[9]*aryubbGL[14];
   aryubbGL[15]=3*aryubbGL[6] + 1./2.*MMt;
   aryubbGL[15]=aryubbGL[8]*aryubbGL[15];
   aryubbGL[16]= - aryubbGL[7]*aryubbGL[8];
   aryubbGL[15]=aryubbGL[15] + 5*aryubbGL[16];
   aryubbGL[15]=aryubbGL[7]*aryubbGL[15];
   aryubbGL[13]= - 103./12. + 32*aryubbGL[13];
   aryubbGL[13]=MMt*aryubbGL[13];
   aryubbGL[13]=aryubbGL[14] + 1./2.*aryubbGL[15] + aryubbGL[13] - 1./2.
   *aryubbGL[6] - 4*aryubbGL[5] + 5*aryubbGL[12] + 4*aryubbGL[11];

      yubbGLret = aryubbGL[13]*aryubbGL[4]*pow(aryubbGL[2],2)*
      aryubbGL[1];
      return yubbGLret;
}
