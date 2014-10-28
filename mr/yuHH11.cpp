#include <HH.hpp>
std::complex<long double>
HH<OS>::my11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHH[19], yuHHret;

    aryuHH[1]=double(nH);
    aryuHH[2]=pow(CW,-1);
    aryuHH[3]=pow(MMH,-1);
    aryuHH[4]=pow(MMZ,-1);
    aryuHH[5]=pow(SW,-1);
    aryuHH[6]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHH[7]=Tsil::A(MMt,mu2);
    aryuHH[8]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHH[9]=Tsil::Aeps(MMt,mu2);
    aryuHH[10]=prottttt0->M(0);
    aryuHH[11]=prottttt0->Vxzuv(0);
    aryuHH[12]=prottttt0->Suxv(0);
   aryuHH[13]= - 2*aryuHH[11] - 3*aryuHH[10];
   aryuHH[14]=pow(aryuHH[5],2);
   aryuHH[15]=aryuHH[14]*aryuHH[13];
   aryuHH[16]=pow(aryuHH[2],2);
   aryuHH[13]=aryuHH[16]*aryuHH[13];
   aryuHH[13]=aryuHH[15] + aryuHH[13];
   aryuHH[15]= - 10 - 7*aryuHH[6];
   aryuHH[15]=aryuHH[6]*aryuHH[15];
   aryuHH[15]=aryuHH[15] - 9 + 4*aryuHH[8];
   aryuHH[17]=aryuHH[14]*aryuHH[15];
   aryuHH[15]=aryuHH[16]*aryuHH[15];
   aryuHH[15]=aryuHH[17] + aryuHH[15];
   aryuHH[15]=aryuHH[3]*aryuHH[15];
   aryuHH[17]=2*aryuHH[11] + aryuHH[10];
   aryuHH[18]=aryuHH[14]*aryuHH[17];
   aryuHH[17]=aryuHH[16]*aryuHH[17];
   aryuHH[17]=aryuHH[18] + aryuHH[17];
   aryuHH[17]=MMt*aryuHH[3]*aryuHH[17];
   aryuHH[13]=8*aryuHH[17] + 2*aryuHH[13] + aryuHH[15];
   aryuHH[13]=MMt*aryuHH[13];
   aryuHH[15]= - aryuHH[6]*aryuHH[7];
   aryuHH[15]=18*aryuHH[15] - aryuHH[12] + 6*aryuHH[9];
   aryuHH[17]=aryuHH[14]*aryuHH[15];
   aryuHH[15]=aryuHH[16]*aryuHH[15];
   aryuHH[15]=aryuHH[17] + aryuHH[15];
   aryuHH[15]=aryuHH[3]*aryuHH[15];
   aryuHH[17]=aryuHH[10]*MMH;
   aryuHH[17]=6 + aryuHH[17];
   aryuHH[18]=4 + 3*aryuHH[6];
   aryuHH[18]=aryuHH[6]*aryuHH[18];
   aryuHH[17]=2*aryuHH[17] + aryuHH[18];
   aryuHH[18]=aryuHH[14]*aryuHH[17];
   aryuHH[17]=aryuHH[16]*aryuHH[17];
   aryuHH[13]=2*aryuHH[13] + 2*aryuHH[15] + aryuHH[18] + aryuHH[17];
   aryuHH[13]=MMt*aryuHH[13];
   aryuHH[15]=pow(aryuHH[7],2);
   aryuHH[17]= - aryuHH[14]*aryuHH[15];
   aryuHH[15]= - aryuHH[16]*aryuHH[15];
   aryuHH[15]=aryuHH[17] + aryuHH[15];
   aryuHH[15]=aryuHH[3]*aryuHH[15];
   aryuHH[17]=aryuHH[6]*aryuHH[7];
   aryuHH[14]=aryuHH[14]*aryuHH[17];
   aryuHH[16]=aryuHH[16]*aryuHH[17];
   aryuHH[14]=4*aryuHH[15] + aryuHH[14] + aryuHH[16];
   aryuHH[13]=6*aryuHH[14] + aryuHH[13];

      yuHHret = 2*aryuHH[13]*aryuHH[4]*aryuHH[1];
      return yuHHret;
}
