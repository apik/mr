#include <HH.hpp>
std::complex<long double>
HH<OS>::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHH[21], yuHHret;

    aryuHH[1]=double(nH);
    aryuHH[2]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHH[3]=pow(CW,-1);
    aryuHH[4]=pow(MMH,-1);
    aryuHH[5]=pow(MMZ,-1);
    aryuHH[6]=pow(SW,-1);
    aryuHH[7]=Tsil::A(MMt,mu2);
    aryuHH[8]=double(boson);
    aryuHH[9]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHH[10]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryuHH[11]=Tsil::B(MMW,MMW,MMH,mu2);
    aryuHH[12]=Tsil::A(MMH,mu2);
    aryuHH[13]=Tsil::A(MMZ,mu2);
    aryuHH[14]=Tsil::A(MMW,mu2);
   aryuHH[15]=1./2.*aryuHH[13];
   aryuHH[16]=9./2.*aryuHH[9] + aryuHH[11];
   aryuHH[16]=MMH*aryuHH[16];
   aryuHH[17]=aryuHH[10]*MMH;
   aryuHH[16]=1./4.*aryuHH[17] + 1./2.*aryuHH[16] + aryuHH[15] + 3./2.*
   aryuHH[12] + aryuHH[14];
   aryuHH[17]=pow(aryuHH[3],2);
   aryuHH[18]=aryuHH[17]*aryuHH[16];
   aryuHH[19]=pow(aryuHH[6],2);
   aryuHH[16]=aryuHH[19]*aryuHH[16];
   aryuHH[16]=aryuHH[18] + aryuHH[16];
   aryuHH[16]=aryuHH[5]*aryuHH[16];
   aryuHH[18]=aryuHH[11]*MMZ;
   aryuHH[20]=aryuHH[10]*MMZ;
   aryuHH[15]=1./2.*aryuHH[20] + aryuHH[18] + aryuHH[14] + aryuHH[15];
   aryuHH[15]=aryuHH[4]*aryuHH[15];
   aryuHH[15]=3*aryuHH[15] - aryuHH[11] - 1./2.*aryuHH[10];
   aryuHH[15]=aryuHH[19]*aryuHH[15];
   aryuHH[18]=aryuHH[13] + aryuHH[20];
   aryuHH[18]=aryuHH[4]*aryuHH[18];
   aryuHH[18]= - aryuHH[10] + 3*aryuHH[18];
   aryuHH[18]=aryuHH[17]*aryuHH[18];
   aryuHH[20]= - aryuHH[4]*aryuHH[11]*MMZ;
   aryuHH[15]=1./2.*aryuHH[16] + aryuHH[15] + 3*aryuHH[20] + 1./2.*
   aryuHH[18];
   aryuHH[15]=aryuHH[8]*aryuHH[15];
   aryuHH[16]=aryuHH[1]*MMt*aryuHH[2];
   aryuHH[18]= - MMt*aryuHH[2];
   aryuHH[18]= - aryuHH[7] + aryuHH[18];
   aryuHH[18]=aryuHH[4]*aryuHH[1]*MMt*aryuHH[18];
   aryuHH[16]=1./2.*aryuHH[16] + 2*aryuHH[18];
   aryuHH[17]=aryuHH[17]*aryuHH[16];
   aryuHH[16]=aryuHH[19]*aryuHH[16];
   aryuHH[16]=aryuHH[17] + aryuHH[16];
   aryuHH[16]=aryuHH[5]*aryuHH[16];

      yuHHret = aryuHH[15] + 3*aryuHH[16];
      return yuHHret;
}
