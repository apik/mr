#include <HH.hpp>
std::complex<long double>
HH<OS>::mygl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHHGL[12], yuHHGLret;

    aryuHHGL[1]=double(boson);
    aryuHHGL[2]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHHGL[3]=pow(SW,-1);
    aryuHHGL[4]=pow(MMW,-1);
    aryuHHGL[5]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHHGL[6]=pow(MMH,-1);
    aryuHHGL[7]=Tsil::A(MMH,mu2);
    aryuHHGL[8]=Tsil::A(MMt,mu2);
    aryuHHGL[9]=std::real(Tsil::B(0,0,MMH,mu2));
   aryuHHGL[10]= - aryuHHGL[6]*aryuHHGL[8];
   aryuHHGL[11]= - MMt*aryuHHGL[5]*aryuHHGL[6];
   aryuHHGL[10]=2*aryuHHGL[11] + 2*aryuHHGL[10] + 1./2.*aryuHHGL[5];
   aryuHHGL[10]=MMt*aryuHHGL[10];
   aryuHHGL[11]=aryuHHGL[9] + 3*aryuHHGL[2];
   aryuHHGL[11]=MMH*aryuHHGL[11];
   aryuHHGL[11]=aryuHHGL[7] + 1./2.*aryuHHGL[11];
   aryuHHGL[10]=1./4.*aryuHHGL[11] + aryuHHGL[10];

      yuHHGLret = 3*aryuHHGL[10]*aryuHHGL[4]*pow(aryuHHGL[3],2)*
      aryuHHGL[1];
      return yuHHGLret;
}
