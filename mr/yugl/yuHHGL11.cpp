#include <HH.hpp>
namespace mr
{
  double HH<OS>::ygl11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryuHHGL[20], yuHHGLret;

    aryuHHGL[1]=double(boson);
    aryuHHGL[2]=pow(SW,-1);
    aryuHHGL[3]=pow(MMH,-1);
    aryuHHGL[4]=pow(MMW,-1);
    aryuHHGL[5]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHHGL[6]=Tsil::A(MMt,mu2);
    aryuHHGL[7]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHHGL[8]=pow(MMt,-1);
    aryuHHGL[9]=Tsil::Aeps(MMt,mu2);
    aryuHHGL[10]=prottttt0->M(0);
    aryuHHGL[11]=prottttt0->Vxzuv(0);
    aryuHHGL[12]=prottttt0->Suxv(0);
    aryuHHGL[13]=2*aryuHHGL[5];
    aryuHHGL[14]=4 + 3*aryuHHGL[5];
    aryuHHGL[13]=aryuHHGL[14]*aryuHHGL[13];
    aryuHHGL[14]=4*aryuHHGL[3];
    aryuHHGL[15]=aryuHHGL[5]*aryuHHGL[6];
    aryuHHGL[16]=6*aryuHHGL[9] - aryuHHGL[12] - 18*aryuHHGL[15];
    aryuHHGL[14]=aryuHHGL[16]*aryuHHGL[14];
    aryuHHGL[16]=4*MMt;
    aryuHHGL[17]=aryuHHGL[16]*aryuHHGL[3];
    aryuHHGL[18]= - 10 - 7*aryuHHGL[5];
    aryuHHGL[18]=aryuHHGL[5]*aryuHHGL[18];
    aryuHHGL[18]= - 25 + aryuHHGL[18];
    aryuHHGL[18]=aryuHHGL[18]*aryuHHGL[17];
    aryuHHGL[19]=pow(Pi,2);
    aryuHHGL[13]=aryuHHGL[18] - 2./3.*aryuHHGL[19] + aryuHHGL[14] + 77./
      2. + aryuHHGL[13];
    aryuHHGL[13]=MMt*aryuHHGL[13];
    aryuHHGL[14]=aryuHHGL[7]*aryuHHGL[3];
    aryuHHGL[18]= - 1 + aryuHHGL[17];
    aryuHHGL[18]=aryuHHGL[11]*aryuHHGL[18];
    aryuHHGL[14]=aryuHHGL[18] + aryuHHGL[14];
    aryuHHGL[18]=pow(MMt,2);
    aryuHHGL[18]=16*aryuHHGL[18];
    aryuHHGL[14]=aryuHHGL[18]*aryuHHGL[14];
    aryuHHGL[18]= - aryuHHGL[6]*aryuHHGL[8];
    aryuHHGL[18]= - 1 + aryuHHGL[18];
    aryuHHGL[18]=aryuHHGL[6]*aryuHHGL[18];
    aryuHHGL[15]=aryuHHGL[18] + 2*aryuHHGL[15];
    aryuHHGL[17]= - 3 + aryuHHGL[17];
    aryuHHGL[17]=MMt*aryuHHGL[17];
    aryuHHGL[17]=MMH + 2*aryuHHGL[17];
    aryuHHGL[16]=aryuHHGL[10]*aryuHHGL[17]*aryuHHGL[16];
    aryuHHGL[13]=aryuHHGL[16] + 6*aryuHHGL[15] + aryuHHGL[13] + 
      aryuHHGL[14];

    yuHHGLret = aryuHHGL[13]*aryuHHGL[4]*pow(aryuHHGL[2],2)*
      aryuHHGL[1];
    return yuHHGLret.real();
  }
} // namespace mr
