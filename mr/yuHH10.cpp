#include <HH.hpp>
std::complex<long double>
HH<OS>::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHH[24], yuHHret;

    aryuHH[1]=double(nH);
    aryuHH[2]=pow(CW,-1);
    aryuHH[3]=pow(MMZ,-1);
    aryuHH[4]=pow(SW,-1);
    aryuHH[5]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHH[6]=pow(MMH,-1);
    aryuHH[7]=Tsil::A(MMt,mu2);
    aryuHH[8]=double(boson);
    aryuHH[9]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHH[10]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryuHH[11]=Tsil::B(MMW,MMW,MMH,mu2);
    aryuHH[12]=Tsil::A(MMZ,mu2);
    aryuHH[13]=Tsil::A(MMW,mu2);
    aryuHH[14]=1/( - MMW + MMH);
    aryuHH[15]=Tsil::A(MMH,mu2);
   aryuHH[16]=9./2.*aryuHH[9] + aryuHH[11] + 1./2.;
   aryuHH[17]=1./4.*MMH;
   aryuHH[16]=aryuHH[16]*aryuHH[17];
   aryuHH[16]=aryuHH[16] + aryuHH[12] + 2*aryuHH[13];
   aryuHH[16]=aryuHH[16]*aryuHH[8];
   aryuHH[17]=aryuHH[10]*aryuHH[8];
   aryuHH[18]=aryuHH[17]*MMH;
   aryuHH[16]=aryuHH[16] + 1./8.*aryuHH[18];
   aryuHH[16]=aryuHH[16]*aryuHH[3];
   aryuHH[18]=2*aryuHH[6];
   aryuHH[18]=aryuHH[18]*MMt;
   aryuHH[18]=aryuHH[18] - 1./2.;
   aryuHH[18]=aryuHH[18]*aryuHH[5];
   aryuHH[18]=aryuHH[18] + 1./4.;
   aryuHH[18]=aryuHH[18]*MMt;
   aryuHH[18]=aryuHH[18] + 1./2.*aryuHH[7];
   aryuHH[19]=aryuHH[1]*aryuHH[3];
   aryuHH[19]=3*aryuHH[19];
   aryuHH[18]=aryuHH[18]*aryuHH[19];
   aryuHH[19]=aryuHH[6]*MMZ;
   aryuHH[20]= - 1 + 3*aryuHH[19];
   aryuHH[17]=1./2.*aryuHH[17];
   aryuHH[17]=aryuHH[20]*aryuHH[17];
   aryuHH[16]=aryuHH[17] + aryuHH[16] - aryuHH[18];
   aryuHH[17]=aryuHH[19] - 1./8.;
   aryuHH[18]= - aryuHH[8]*aryuHH[17];
   aryuHH[18]=aryuHH[18] + aryuHH[16];
   aryuHH[18]=aryuHH[18]*pow(aryuHH[2],2);
   aryuHH[21]=pow(aryuHH[4],2);
   aryuHH[22]=aryuHH[13] - aryuHH[12];
   aryuHH[22]=aryuHH[3]*aryuHH[22]*aryuHH[21];
   aryuHH[23]= - aryuHH[13] + aryuHH[15];
   aryuHH[23]=aryuHH[14]*aryuHH[23];
   aryuHH[22]=aryuHH[23] + aryuHH[22];
   aryuHH[20]=aryuHH[11]*aryuHH[20];
   aryuHH[17]= - 3*aryuHH[17] + aryuHH[20] + 3./4.*aryuHH[22];
   aryuHH[17]=aryuHH[8]*aryuHH[17];
   aryuHH[16]=aryuHH[17] + aryuHH[16];
   aryuHH[16]=aryuHH[16]*aryuHH[21];
   aryuHH[17]=2 - 3*aryuHH[11];
   aryuHH[17]=aryuHH[8]*aryuHH[19]*aryuHH[17];

      yuHHret = aryuHH[16] + aryuHH[17] + aryuHH[18];
      return yuHHret;
}
