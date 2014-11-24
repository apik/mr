#include <bb.hpp>
std::complex<long double>
bb::y10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubb[28], yubbret;

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
   aryubb[17]=1./2.*MMt;
   aryubb[18]=1./2.*MMb;
   aryubb[19]=aryubb[17] + aryubb[6] - aryubb[18];
   aryubb[20]=3*MMb;
   aryubb[20]=aryubb[20]*aryubb[15]*aryubb[1];
   aryubb[19]=aryubb[19]*aryubb[20];
   aryubb[21]=MMt*aryubb[1];
   aryubb[21]=aryubb[21] - aryubb[12];
   aryubb[22]= - aryubb[11] + aryubb[6] + 7*aryubb[13];
   aryubb[20]=aryubb[20] + 1;
   aryubb[20]=aryubb[20]*aryubb[14];
   aryubb[23]=1./4.*MMH;
   aryubb[24]=aryubb[1]*aryubb[6];
   aryubb[19]= - 3*aryubb[24] + aryubb[20] + aryubb[23] - aryubb[19] + 
   1./2.*aryubb[22] - 3./2.*aryubb[21];
   aryubb[18]=aryubb[18] - MMt;
   aryubb[20]=aryubb[10]*aryubb[18];
   aryubb[20]=aryubb[20] + aryubb[19];
   aryubb[20]=aryubb[4]*aryubb[20];
   aryubb[21]=aryubb[13] - aryubb[6];
   aryubb[22]= - MMt*aryubb[21];
   aryubb[24]=pow(MMt,2);
   aryubb[25]=aryubb[10]*aryubb[24];
   aryubb[22]=aryubb[25] + aryubb[22];
   aryubb[22]=aryubb[4]*aryubb[22];
   aryubb[25]=aryubb[8]*MMZ;
   aryubb[25]=aryubb[25] + aryubb[12];
   aryubb[26]=aryubb[25] - aryubb[14];
   aryubb[22]= - 5./9.*aryubb[26] + aryubb[22];
   aryubb[22]=aryubb[9]*aryubb[22];
   aryubb[27]= - 1./2. + aryubb[8];
   aryubb[20]=1./2.*aryubb[22] + 17./18.*aryubb[27] + aryubb[20];
   aryubb[22]=aryubb[23] - MMb;
   aryubb[22]=aryubb[22]*aryubb[7];
   aryubb[23]= - aryubb[4]*aryubb[22];
   aryubb[20]=1./2.*aryubb[20] + aryubb[23];
   aryubb[20]=aryubb[20]*pow(aryubb[3],2);
   aryubb[23]=pow(aryubb[5],2);
   aryubb[27]=aryubb[23]*aryubb[4];
   aryubb[22]=aryubb[27]*aryubb[22];
   aryubb[20]=aryubb[22] - aryubb[20];
   aryubb[19]=aryubb[4]*aryubb[19];
   aryubb[22]=aryubb[13] - aryubb[12];
   aryubb[22]=aryubb[22]*aryubb[27];
   aryubb[27]= - aryubb[13] + aryubb[11];
   aryubb[27]=aryubb[16]*aryubb[27];
   aryubb[27]= - 1./2. + aryubb[27];
   aryubb[27]=3*aryubb[27] + aryubb[8];
   aryubb[18]=aryubb[4]*aryubb[18];
   aryubb[18]=1./2. + aryubb[18];
   aryubb[18]=aryubb[10]*aryubb[18];
   aryubb[18]=aryubb[18] + 3./2.*aryubb[22] + 1./2.*aryubb[27] + 
   aryubb[19];
   aryubb[19]=1./4.*aryubb[23];
   aryubb[18]=aryubb[19]*aryubb[18];
   aryubb[22]= - aryubb[4]*aryubb[17];
   aryubb[22]=aryubb[22] - 1;
   aryubb[21]=aryubb[21]*aryubb[22];
   aryubb[21]=aryubb[21] - 1./2.*aryubb[26];
   aryubb[19]=aryubb[21]*aryubb[19];
   aryubb[21]=aryubb[4]*aryubb[24];
   aryubb[17]=1./2.*aryubb[21] - MMZ + aryubb[17];
   aryubb[17]=aryubb[17]*aryubb[23];
   aryubb[17]=MMZ + aryubb[17];
   aryubb[17]=aryubb[10]*aryubb[17];
   aryubb[21]=2*aryubb[14] + aryubb[25];
   aryubb[17]=1./4.*aryubb[17] + 1./9.*aryubb[21] + aryubb[19];
   aryubb[17]=aryubb[9]*aryubb[17];
   aryubb[19]= - 1 + aryubb[8];
   aryubb[17]=aryubb[17] + 2./9.*aryubb[19] + aryubb[18] - 1./2.*
   aryubb[20];

      yubbret = aryubb[17]*aryubb[2];
      return yubbret;
}
