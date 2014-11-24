#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::y10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuZZ[30], yuZZret;

    aryuZZ[1]=double(nH);
    aryuZZ[2]=double(boson);
    aryuZZ[3]=pow(CW,-1);
    aryuZZ[4]=pow(MMZ,-1);
    aryuZZ[5]=pow(SW,-1);
    aryuZZ[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuZZ[7]=Tsil::B(MMb,MMb,MMZ,mu2);
    aryuZZ[8]=Tsil::B(0,0,MMZ,mu2);
    aryuZZ[9]=Tsil::A(MMt,mu2);
    aryuZZ[10]=Tsil::A(MMb,mu2);
    aryuZZ[11]=double(nL + nH);
    aryuZZ[12]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryuZZ[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryuZZ[14]=Tsil::A(MMH,mu2);
    aryuZZ[15]=Tsil::A(MMZ,mu2);
    aryuZZ[16]=Tsil::A(MMW,mu2);
    aryuZZ[17]=1/( - MMb + MMt);
    aryuZZ[18]=1/( - MMW + MMH);
   aryuZZ[19]=pow(aryuZZ[3],2);
   aryuZZ[20]=pow(aryuZZ[5],2);
   aryuZZ[21]=5./9.*aryuZZ[19] + aryuZZ[20] - 8./9.;
   aryuZZ[22]=aryuZZ[19] + aryuZZ[20];
   aryuZZ[23]=aryuZZ[22]*aryuZZ[17];
   aryuZZ[24]=MMb*aryuZZ[23];
   aryuZZ[25]= - 3./4.*MMt - 3./2.*aryuZZ[9];
   aryuZZ[23]=aryuZZ[25]*aryuZZ[23];
   aryuZZ[25]=1./2.*aryuZZ[20];
   aryuZZ[26]=aryuZZ[25] + 17./18.*aryuZZ[19];
   aryuZZ[27]= - 8./9. - aryuZZ[26];
   aryuZZ[27]=aryuZZ[7]*aryuZZ[27];
   aryuZZ[23]=3./4.*aryuZZ[24] + aryuZZ[27] + aryuZZ[23] + aryuZZ[21];
   aryuZZ[23]=MMb*aryuZZ[23];
   aryuZZ[21]=aryuZZ[21] + 3./2.*aryuZZ[24];
   aryuZZ[21]=aryuZZ[10]*aryuZZ[21];
   aryuZZ[21]=aryuZZ[21] + aryuZZ[23];
   aryuZZ[21]=aryuZZ[21]*aryuZZ[1]*aryuZZ[4];
   aryuZZ[23]= - aryuZZ[25] - 32./9. + 7./18.*aryuZZ[19];
   aryuZZ[24]=aryuZZ[6]*aryuZZ[23];
   aryuZZ[24]=aryuZZ[24] + 41./36.*aryuZZ[19] - 32./9. + 1./4.*
   aryuZZ[20];
   aryuZZ[24]=MMt*aryuZZ[24];
   aryuZZ[23]=aryuZZ[9]*aryuZZ[23];
   aryuZZ[23]=aryuZZ[23] + aryuZZ[24];
   aryuZZ[23]=aryuZZ[4]*aryuZZ[23];
   aryuZZ[24]= - 11./9.*aryuZZ[19] + 20./9. - aryuZZ[20];
   aryuZZ[24]=aryuZZ[8]*aryuZZ[24];
   aryuZZ[25]=5./18.*aryuZZ[19] - 4./9. + aryuZZ[25];
   aryuZZ[25]=aryuZZ[7]*aryuZZ[25];
   aryuZZ[26]= - 16./9. + aryuZZ[26];
   aryuZZ[26]=aryuZZ[6]*aryuZZ[26];
   aryuZZ[23]=aryuZZ[23] + aryuZZ[26] + aryuZZ[25] + aryuZZ[24];
   aryuZZ[23]=aryuZZ[1]*aryuZZ[23];
   aryuZZ[24]=aryuZZ[22]*aryuZZ[14];
   aryuZZ[25]= - aryuZZ[15]*aryuZZ[22];
   aryuZZ[22]=aryuZZ[22]*MMH;
   aryuZZ[26]=aryuZZ[12]*aryuZZ[22];
   aryuZZ[25]=aryuZZ[24] + aryuZZ[25] + aryuZZ[26];
   aryuZZ[25]=aryuZZ[4]*MMH*aryuZZ[25];
   aryuZZ[26]=11./3. - 3*aryuZZ[20];
   aryuZZ[26]=aryuZZ[26]*aryuZZ[20];
   aryuZZ[26]=aryuZZ[26] + 11./3.*aryuZZ[19];
   aryuZZ[26]=aryuZZ[15]*aryuZZ[26];
   aryuZZ[24]=aryuZZ[24] - aryuZZ[26];
   aryuZZ[26]=aryuZZ[12] + 1./8.;
   aryuZZ[22]= - aryuZZ[26]*aryuZZ[22];
   aryuZZ[22]=1./12.*aryuZZ[25] + 1./3.*aryuZZ[22] - 1./4.*aryuZZ[24];
   aryuZZ[22]=aryuZZ[4]*aryuZZ[22];
   aryuZZ[24]=pow(CW,2);
   aryuZZ[24]=4*aryuZZ[24];
   aryuZZ[25]=aryuZZ[24] + 1./12.*aryuZZ[19] + 29./3. - 33./4.*
   aryuZZ[20];
   aryuZZ[25]=aryuZZ[13]*aryuZZ[25];
   aryuZZ[26]=3./4.*aryuZZ[18];
   aryuZZ[27]= - aryuZZ[20]*aryuZZ[26];
   aryuZZ[28]= - 1 + 3./4.*aryuZZ[20];
   aryuZZ[28]=aryuZZ[28]*aryuZZ[20];
   aryuZZ[29]=5./3.*aryuZZ[19];
   aryuZZ[28]=aryuZZ[29] + 4 + aryuZZ[28];
   aryuZZ[28]=aryuZZ[4]*aryuZZ[28];
   aryuZZ[27]=aryuZZ[27] + aryuZZ[28];
   aryuZZ[27]=aryuZZ[16]*aryuZZ[27];
   aryuZZ[26]=aryuZZ[14]*aryuZZ[26];
   aryuZZ[26]=aryuZZ[26] + aryuZZ[12];
   aryuZZ[26]=aryuZZ[20]*aryuZZ[26];
   aryuZZ[28]=5./72. + aryuZZ[12];
   aryuZZ[19]=aryuZZ[28]*aryuZZ[19];
   aryuZZ[28]=aryuZZ[29] + aryuZZ[20] - 8./3.;
   aryuZZ[29]=aryuZZ[8] - 1./3.;
   aryuZZ[28]=aryuZZ[11]*aryuZZ[29]*aryuZZ[28];
   aryuZZ[20]=8 - 209./24.*aryuZZ[20];
   aryuZZ[19]=4./3.*aryuZZ[28] + aryuZZ[27] + aryuZZ[23] + aryuZZ[22]
    + aryuZZ[25] + aryuZZ[24] + aryuZZ[19] + 1./3.*aryuZZ[20] + 
   aryuZZ[26] + aryuZZ[21];

      yuZZret = aryuZZ[19]*aryuZZ[2];
      return yuZZret;
}
