#include <HH.hpp>
std::complex<long double>
HH<OS>::mygl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHHGL[16], yuHHGLret;

    aryuHHGL[1]=double(boson);
    aryuHHGL[2]=pow(SW,-1);
    aryuHHGL[3]=pow(MMH,-1);
    aryuHHGL[4]=pow(MMW,-1);
    aryuHHGL[5]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHHGL[6]=Tsil::A(MMt,mu2);
    aryuHHGL[7]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHHGL[8]=Tsil::Aeps(MMt,mu2);
    aryuHHGL[9]=prottttt0->M(0);
    aryuHHGL[10]=prottttt0->Vxzuv(0);
    aryuHHGL[11]=prottttt0->Suxv(0);
   aryuHHGL[12]= - 2*aryuHHGL[10] - 3*aryuHHGL[9];
   aryuHHGL[13]= - 10 - 7*aryuHHGL[5];
   aryuHHGL[13]=aryuHHGL[5]*aryuHHGL[13];
   aryuHHGL[13]=aryuHHGL[13] - 9 + 4*aryuHHGL[7];
   aryuHHGL[13]=aryuHHGL[3]*aryuHHGL[13];
   aryuHHGL[14]=2*aryuHHGL[10] + aryuHHGL[9];
   aryuHHGL[14]=MMt*aryuHHGL[3]*aryuHHGL[14];
   aryuHHGL[12]=8*aryuHHGL[14] + 2*aryuHHGL[12] + aryuHHGL[13];
   aryuHHGL[12]=MMt*aryuHHGL[12];
   aryuHHGL[13]=aryuHHGL[9]*MMH;
   aryuHHGL[13]=6 + aryuHHGL[13];
   aryuHHGL[14]=4 + 3*aryuHHGL[5];
   aryuHHGL[14]=aryuHHGL[5]*aryuHHGL[14];
   aryuHHGL[15]= - aryuHHGL[5]*aryuHHGL[6];
   aryuHHGL[15]=18*aryuHHGL[15] - aryuHHGL[11] + 6*aryuHHGL[8];
   aryuHHGL[15]=aryuHHGL[3]*aryuHHGL[15];
   aryuHHGL[12]=2*aryuHHGL[12] + 2*aryuHHGL[15] + 2*aryuHHGL[13] + 
   aryuHHGL[14];
   aryuHHGL[12]=MMt*aryuHHGL[12];
   aryuHHGL[13]=aryuHHGL[5]*aryuHHGL[6];
   aryuHHGL[14]= - aryuHHGL[3]*pow(aryuHHGL[6],2);
   aryuHHGL[13]=aryuHHGL[13] + 4*aryuHHGL[14];
   aryuHHGL[12]=6*aryuHHGL[13] + aryuHHGL[12];

      yuHHGLret = 2*aryuHHGL[12]*aryuHHGL[4]*pow(aryuHHGL[2],2)*
      aryuHHGL[1];
      return yuHHGLret;
}
