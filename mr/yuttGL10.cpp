#include <tt.hpp>
std::complex<long double>
tt::mygl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[11], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[3]=pow(SW,-1);
    aryuttGL[4]=pow(MMW,-1);
    aryuttGL[5]=Tsil::A(MMH,mu2);
    aryuttGL[6]=Tsil::A(MMt,mu2);
    aryuttGL[7]=pow(MMH,-1);
    aryuttGL[8]=std::real(Tsil::B(0,0,MMt,mu2));
   aryuttGL[9]= - aryuttGL[2]*MMH;
   aryuttGL[9]=1./2.*aryuttGL[9] + aryuttGL[5] + aryuttGL[6];
   aryuttGL[10]= - aryuttGL[6]*aryuttGL[7];
   aryuttGL[10]=1./2.*aryuttGL[2] + 1./8.*aryuttGL[8] + 3*aryuttGL[10];
   aryuttGL[10]=MMt*aryuttGL[10];
   aryuttGL[9]=1./4.*aryuttGL[9] + aryuttGL[10];

      yuttGLret = aryuttGL[9]*aryuttGL[4]*pow(aryuttGL[3],2)*
      aryuttGL[1];
      return yuttGLret;
}
