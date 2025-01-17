#include <HH.hpp>
namespace mr
{
  double HH<OS>::ygl10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryuHHGL[11], yuHHGLret;

    aryuHHGL[1]=double(boson);
    aryuHHGL[2]=pow(SW,-1);
    aryuHHGL[3]=pow(MMW,-1);
    aryuHHGL[4]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHHGL[5]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHHGL[6]=pow(MMH,-1);
    aryuHHGL[7]=Tsil::A(MMt,mu2);
    aryuHHGL[8]=std::real(Tsil::B(0,0,MMH,mu2));
    aryuHHGL[9]=MMt*aryuHHGL[6];
    aryuHHGL[9]=1./2. - 2*aryuHHGL[9];
    aryuHHGL[9]=aryuHHGL[5]*aryuHHGL[9];
    aryuHHGL[9]= - 1./4. + aryuHHGL[9];
    aryuHHGL[9]=MMt*aryuHHGL[9];
    aryuHHGL[9]=aryuHHGL[9] - 1./2.*aryuHHGL[7];
    aryuHHGL[10]=3*aryuHHGL[8] + 1 + 9*aryuHHGL[4];
    aryuHHGL[10]=MMH*aryuHHGL[10];
    aryuHHGL[9]=3*aryuHHGL[9] + 1./8.*aryuHHGL[10];

    yuHHGLret = aryuHHGL[9]*aryuHHGL[3]*pow(aryuHHGL[2],2)*
      aryuHHGL[1];
    return yuHHGLret.real();
  }
} // namespace mr
