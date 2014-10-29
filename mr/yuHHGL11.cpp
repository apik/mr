#include <HH.hpp>
std::complex<long double>
HH<OS>::mygl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHHGL[19], yuHHGLret;

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
   aryuHHGL[13]=MMt*aryuHHGL[3];
   aryuHHGL[14]=2*aryuHHGL[11];
   aryuHHGL[15]=aryuHHGL[14] + aryuHHGL[10];
   aryuHHGL[15]=aryuHHGL[15]*aryuHHGL[13];
   aryuHHGL[14]= - aryuHHGL[14] - 3*aryuHHGL[10];
   aryuHHGL[16]= - 25 + 4*aryuHHGL[7];
   aryuHHGL[16]=aryuHHGL[3]*aryuHHGL[16];
   aryuHHGL[17]= - 10 - 7*aryuHHGL[5];
   aryuHHGL[17]=aryuHHGL[5]*aryuHHGL[3]*aryuHHGL[17];
   aryuHHGL[14]=8*aryuHHGL[15] + aryuHHGL[17] + 2*aryuHHGL[14] + 
   aryuHHGL[16];
   aryuHHGL[15]=4*MMt;
   aryuHHGL[14]=aryuHHGL[14]*aryuHHGL[15];
   aryuHHGL[15]=6*aryuHHGL[9] - aryuHHGL[12];
   aryuHHGL[15]=aryuHHGL[3]*aryuHHGL[15];
   aryuHHGL[16]=MMH*aryuHHGL[10];
   aryuHHGL[15]=aryuHHGL[15] + aryuHHGL[16];
   aryuHHGL[16]=2*aryuHHGL[5];
   aryuHHGL[17]=4 + 3*aryuHHGL[5];
   aryuHHGL[17]=aryuHHGL[17]*aryuHHGL[16];
   aryuHHGL[18]=pow(Pi,2);
   aryuHHGL[14]=aryuHHGL[14] + aryuHHGL[17] - 2./3.*aryuHHGL[18] + 77./
   2. + 4*aryuHHGL[15];
   aryuHHGL[14]=aryuHHGL[14]*MMt;
   aryuHHGL[13]=aryuHHGL[5]*aryuHHGL[13];
   aryuHHGL[15]= - aryuHHGL[6]*aryuHHGL[8];
   aryuHHGL[13]=aryuHHGL[15] - 12*aryuHHGL[13] - 1 + aryuHHGL[16];
   aryuHHGL[13]=aryuHHGL[6]*aryuHHGL[13];
   aryuHHGL[13]=6*aryuHHGL[13] + aryuHHGL[14];

      yuHHGLret = aryuHHGL[13]*aryuHHGL[4]*pow(aryuHHGL[2],2)*
      aryuHHGL[1];
      return yuHHGLret;
}
