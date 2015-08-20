#include <bb.hpp>
namespace mr
{
  long double bb<OS>::ygl11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> aryubbGL[14], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMW,-1);
    aryubbGL[4]=Tsil::I2(0,0,MMt,mu2);
    aryubbGL[5]=Tsil::A(MMt,mu2);
    aryubbGL[6]=pow(MMt,-1);
    aryubbGL[7]=Tsil::A(MMb,mu2);
    aryubbGL[8]=pow(MMb,-1);
    aryubbGL[9]=Tsil::Aeps(MMt,mu2);
    aryubbGL[10]=Tsil::Aeps(MMb,mu2);
    aryubbGL[11]=aryubbGL[5]*aryubbGL[6];
    aryubbGL[11]=9 - 7*aryubbGL[11];
    aryubbGL[11]=aryubbGL[5]*aryubbGL[11];
    aryubbGL[12]=1./2.*MMH - 5./2.*MMt - 3*aryubbGL[5] - 5*aryubbGL[7];
    aryubbGL[12]=aryubbGL[8]*aryubbGL[7]*aryubbGL[12];
    aryubbGL[11]=aryubbGL[11] + aryubbGL[12];
    aryubbGL[12]=aryubbGL[9] - aryubbGL[4];
    aryubbGL[13]=pow(Pi,2);
    aryubbGL[13]= - 5./2. - aryubbGL[13];
    aryubbGL[13]=MMt*aryubbGL[13];
    aryubbGL[11]= - 1./12.*MMH + 1./3.*aryubbGL[13] + 5*aryubbGL[10] + 1.
      /2.*aryubbGL[11] + 4*aryubbGL[12];

    yubbGLret = aryubbGL[11]*aryubbGL[3]*pow(aryubbGL[2],2)*
      aryubbGL[1];
    return yubbGLret.real();
  }
} // namespace mr
