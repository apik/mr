#include <tt.hpp>
namespace mr
{
  long double tt<OS>::ygl10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> aryuttGL[10], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=pow(SW,-1);
    aryuttGL[3]=pow(MMW,-1);
    aryuttGL[4]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[5]=Tsil::A(MMH,mu2);
    aryuttGL[6]=Tsil::A(MMt,mu2);
    aryuttGL[7]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuttGL[8]=1./2. - aryuttGL[4];
    aryuttGL[8]=MMH*aryuttGL[8];
    aryuttGL[8]=aryuttGL[8] - aryuttGL[5];
    aryuttGL[9]= - 3 + aryuttGL[7];
    aryuttGL[9]=1./4.*aryuttGL[9] + aryuttGL[4];
    aryuttGL[9]=MMt*aryuttGL[9];
    aryuttGL[8]=aryuttGL[9] - aryuttGL[6] + 1./4.*aryuttGL[8];

    yuttGLret = 1./2.*aryuttGL[8]*aryuttGL[3]*pow(aryuttGL[2],2)*
      aryuttGL[1];
    return yuttGLret.real();
  }
} // namespace mr
