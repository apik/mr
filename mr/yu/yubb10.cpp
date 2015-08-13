#include <bb.hpp>
long double bb<OS>::y10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubb[27], yubbret;

    aryubb[1]=double(nH);
    aryubb[2]=double(boson);
    aryubb[3]=pow(CW,-1);
    aryubb[4]=pow(MMZ,-1);
    aryubb[5]=pow(SW,-1);
    aryubb[6]=Tsil::A(MMt,mu2);
    aryubb[7]=Tsil::B(MMH,MMb,MMb,mu2);
    aryubb[8]=Tsil::B(MMZ,MMb,MMb,mu2);
    aryubb[9]=pow(MMb,-1);
    aryubb[10]=Tsil::B(MMW,MMt,MMb,mu2);
    aryubb[11]=Tsil::A(MMH,mu2);
    aryubb[12]=Tsil::A(MMZ,mu2);
    aryubb[13]=Tsil::A(MMW,mu2);
    aryubb[14]=Tsil::A(MMb,mu2);
    aryubb[15]=1/( - MMb + MMt);
    aryubb[16]=1/( - MMW + MMH);
   aryubb[17]=aryubb[11] - aryubb[13];
   aryubb[18]=3./2.*aryubb[16];
   aryubb[17]=aryubb[18]*aryubb[17];
   aryubb[18]=MMZ*aryubb[9];
   aryubb[19]=aryubb[18]*aryubb[10];
   aryubb[20]=1./2.*MMt;
   aryubb[21]=aryubb[20]*aryubb[9]*aryubb[10];
   aryubb[22]= - 3./2. + aryubb[10];
   aryubb[23]=aryubb[13] - aryubb[6];
   aryubb[24]= - aryubb[9]*aryubb[23];
   aryubb[25]=1./2.*aryubb[9];
   aryubb[26]= - aryubb[12]*aryubb[25];
   aryubb[17]=aryubb[21] + aryubb[26] - aryubb[19] + 1./2.*aryubb[22]
    + aryubb[24] + aryubb[17];
   aryubb[22]=pow(aryubb[5],2);
   aryubb[17]=aryubb[17]*aryubb[22];
   aryubb[17]=aryubb[19] + aryubb[17];
   aryubb[19]=aryubb[23]*aryubb[25];
   aryubb[19]= - aryubb[19] + aryubb[21] - aryubb[10];
   aryubb[19]=aryubb[19]*aryubb[20];
   aryubb[21]=MMb*aryubb[10];
   aryubb[21]= - aryubb[11] + aryubb[6] + 1./2.*MMH + aryubb[21] + 7*
   aryubb[13];
   aryubb[23]= - MMb + 1./4.*MMH;
   aryubb[23]=aryubb[23]*aryubb[7];
   aryubb[19]=3./4.*aryubb[12] + aryubb[19] - aryubb[23] + 1./4.*
   aryubb[21];
   aryubb[21]=pow(aryubb[3],2);
   aryubb[23]=aryubb[19]*aryubb[21];
   aryubb[24]=aryubb[13] - aryubb[12];
   aryubb[24]=aryubb[24]*aryubb[22];
   aryubb[19]=3./4.*aryubb[24] + aryubb[19];
   aryubb[19]=aryubb[19]*aryubb[22];
   aryubb[24]=aryubb[15]*MMb;
   aryubb[25]=aryubb[24] + 1;
   aryubb[20]=aryubb[25]*aryubb[20];
   aryubb[25]= - aryubb[6] + 1./2.*MMb;
   aryubb[25]=aryubb[25]*aryubb[24];
   aryubb[20]=aryubb[6] + aryubb[20] - aryubb[25];
   aryubb[25]=aryubb[21] + aryubb[22];
   aryubb[20]= - aryubb[20]*aryubb[25];
   aryubb[25]=aryubb[25]*aryubb[14];
   aryubb[24]=aryubb[24]*aryubb[25];
   aryubb[20]=aryubb[24] + aryubb[20];
   aryubb[20]=aryubb[1]*aryubb[20];
   aryubb[19]=3./2.*aryubb[20] + 1./2.*aryubb[25] + aryubb[23] + 
   aryubb[19];
   aryubb[19]=aryubb[4]*aryubb[19];
   aryubb[20]=17 - 5*aryubb[18];
   aryubb[20]=aryubb[20]*aryubb[21];
   aryubb[20]=1./8.*aryubb[20] + 2 + aryubb[18];
   aryubb[22]=1./8.*aryubb[22];
   aryubb[18]=1 - aryubb[18];
   aryubb[18]=aryubb[18]*aryubb[22];
   aryubb[18]=1./9.*aryubb[20] + aryubb[18];
   aryubb[18]=aryubb[8]*aryubb[18];
   aryubb[20]=aryubb[12]*aryubb[9];
   aryubb[23]= - 17./2. - 5*aryubb[20];
   aryubb[23]=aryubb[23]*aryubb[21];
   aryubb[21]=2 + 5./8.*aryubb[21];
   aryubb[21]=1./9.*aryubb[21] + aryubb[22];
   aryubb[21]=aryubb[14]*aryubb[9]*aryubb[21];
   aryubb[17]=aryubb[18] + 1./2.*aryubb[19] + aryubb[21] + 1./72.*
   aryubb[23] + 1./9.*aryubb[20] - 2./9. + 1./4.*aryubb[17];

      yubbret = aryubb[17]*aryubb[2];
      return yubbret.real();
}
