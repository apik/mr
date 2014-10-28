#include <bb.hpp>
std::complex<long double>
bb::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubb[20], yubbret;

    aryubb[1]=double(boson);
    aryubb[2]=pow(CW,-1);
    aryubb[3]=pow(MMH,-1);
    aryubb[4]=pow(MMZ,-1);
    aryubb[5]=pow(SW,-1);
    aryubb[6]=Tsil::A(MMH,mu2);
    aryubb[7]=Tsil::A(MMZ,mu2);
    aryubb[8]=Tsil::A(MMW,mu2);
    aryubb[9]=Tsil::A(MMt,mu2);
    aryubb[10]=Tsil::A(MMb,mu2);
    aryubb[11]=pow(MMb,-1);
    aryubb[12]=1/(MMt - MMW);
   aryubb[13]=aryubb[12]*aryubb[8];
   aryubb[13]= - 1./2. + aryubb[13];
   aryubb[14]=aryubb[8] + 1./2.*aryubb[7];
   aryubb[14]=aryubb[3]*aryubb[14];
   aryubb[15]=aryubb[8] - aryubb[9];
   aryubb[15]=aryubb[12]*aryubb[15];
   aryubb[15]=1 + aryubb[15];
   aryubb[15]=aryubb[12]*aryubb[15];
   aryubb[15]=1./4.*aryubb[15] + aryubb[3];
   aryubb[15]=MMZ*aryubb[15];
   aryubb[13]=aryubb[15] + 1./4.*aryubb[13] + aryubb[14];
   aryubb[14]=3*aryubb[6] + 1./2.*MMt;
   aryubb[15]=aryubb[14] + 3*aryubb[9];
   aryubb[16]= - 3*aryubb[3]*aryubb[9]*MMt;
   aryubb[15]=1./8.*aryubb[15] + aryubb[16];
   aryubb[15]=aryubb[4]*aryubb[15];
   aryubb[13]=3./2.*aryubb[13] + aryubb[15];
   aryubb[13]=aryubb[13]*pow(aryubb[5],2);
   aryubb[15]=aryubb[3]*aryubb[7];
   aryubb[15]= - 31./36. + 3*aryubb[15];
   aryubb[17]=pow(aryubb[2],2);
   aryubb[15]=aryubb[17]*aryubb[15];
   aryubb[18]= - aryubb[8] + aryubb[9];
   aryubb[18]=aryubb[12]*aryubb[18];
   aryubb[18]= - 1 + aryubb[18];
   aryubb[18]=aryubb[12]*aryubb[18];
   aryubb[19]=aryubb[17]*aryubb[3];
   aryubb[18]=1./2.*aryubb[19] + 3./8.*aryubb[18] - aryubb[3];
   aryubb[18]=MMZ*aryubb[18];
   aryubb[19]= - 1./3.*aryubb[7];
   aryubb[14]=3./4.*aryubb[9] + 1./4.*aryubb[14] + aryubb[19];
   aryubb[14]=1./2.*aryubb[14] + aryubb[16];
   aryubb[14]=aryubb[17]*aryubb[14];
   aryubb[14]=aryubb[19] + aryubb[14];
   aryubb[14]=aryubb[4]*aryubb[14];
   aryubb[16]=aryubb[10]*aryubb[11];
   aryubb[16]= - 1./2. + aryubb[16];
   aryubb[13]=aryubb[13] + aryubb[14] + aryubb[18] + 1./3.*aryubb[16]
    + 1./4.*aryubb[15];

      yubbret = aryubb[13]*aryubb[1];
      return yubbret;
}
