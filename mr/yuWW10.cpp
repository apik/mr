#include <WW.hpp>
std::complex<long double>
WW<OS>::y10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuWW[27], yuWWret;

    aryuWW[1]=double(nH);
    aryuWW[2]=double(boson);
    aryuWW[3]=pow(CW,-1);
    aryuWW[4]=pow(MMZ,-1);
    aryuWW[5]=pow(SW,-1);
    aryuWW[6]=Tsil::B(MMt,MMb,MMW,mu2);
    aryuWW[7]=Tsil::B(0,0,MMW,mu2);
    aryuWW[8]=Tsil::A(MMt,mu2);
    aryuWW[9]=Tsil::A(MMb,mu2);
    aryuWW[10]=double(nL + nH);
    aryuWW[11]=Tsil::B(MMW,MMH,MMW,mu2);
    aryuWW[12]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryuWW[13]=Tsil::A(MMH,mu2);
    aryuWW[14]=Tsil::A(MMZ,mu2);
    aryuWW[15]=Tsil::A(MMW,mu2);
    aryuWW[16]=1/( - MMb + MMt);
    aryuWW[17]=1/( - MMW + MMH);
   aryuWW[18]=aryuWW[13] - aryuWW[15];
   aryuWW[19]=pow(aryuWW[5],2);
   aryuWW[20]=aryuWW[18]*aryuWW[19];
   aryuWW[21]=pow(aryuWW[3],2);
   aryuWW[22]=aryuWW[21] + 1;
   aryuWW[22]=aryuWW[22]*aryuWW[21];
   aryuWW[18]=aryuWW[18]*aryuWW[22];
   aryuWW[22]=aryuWW[22] + aryuWW[19];
   aryuWW[23]=MMH*aryuWW[11]*aryuWW[22];
   aryuWW[18]=aryuWW[23] + aryuWW[20] + aryuWW[18];
   aryuWW[18]=MMH*aryuWW[18];
   aryuWW[23]=1./2.*MMt;
   aryuWW[24]=aryuWW[23] - MMb;
   aryuWW[24]=aryuWW[24]*MMt;
   aryuWW[25]=pow(MMb,2);
   aryuWW[24]=aryuWW[24] + 1./2.*aryuWW[25];
   aryuWW[24]=aryuWW[24]*aryuWW[6];
   aryuWW[25]=aryuWW[9] - aryuWW[8];
   aryuWW[26]=MMt - MMb;
   aryuWW[26]=aryuWW[26]*aryuWW[25];
   aryuWW[24]=aryuWW[24] - 1./2.*aryuWW[26];
   aryuWW[22]= - aryuWW[1]*aryuWW[24]*aryuWW[22];
   aryuWW[18]=1./12.*aryuWW[18] + aryuWW[22];
   aryuWW[18]=aryuWW[4]*aryuWW[18];
   aryuWW[22]= - aryuWW[23] + 1./2.*MMb + aryuWW[25];
   aryuWW[23]=3./2.*aryuWW[16];
   aryuWW[22]=aryuWW[23]*aryuWW[22];
   aryuWW[22]=aryuWW[22] + 1;
   aryuWW[22]=MMb*aryuWW[22];
   aryuWW[23]=MMt + MMb;
   aryuWW[24]=1./2.*aryuWW[6];
   aryuWW[23]=aryuWW[23]*aryuWW[24];
   aryuWW[22]= - aryuWW[23] + aryuWW[9] + 1./4.*MMt - 1./2.*aryuWW[8]
    + aryuWW[22];
   aryuWW[22]=aryuWW[1]*aryuWW[22];
   aryuWW[23]=aryuWW[11] + 1./8.;
   aryuWW[23]= - MMH*aryuWW[23];
   aryuWW[22]=aryuWW[22] + 1./3.*aryuWW[23];
   aryuWW[23]=aryuWW[19] + aryuWW[21];
   aryuWW[22]=aryuWW[23]*aryuWW[22];
   aryuWW[23]=aryuWW[14] - aryuWW[15];
   aryuWW[24]=aryuWW[23]*aryuWW[19];
   aryuWW[25]=1./2.*aryuWW[13];
   aryuWW[24]= - 3./2.*aryuWW[24] - aryuWW[25] + 5./6.*aryuWW[15] - 
   aryuWW[14];
   aryuWW[24]=aryuWW[24]*aryuWW[19];
   aryuWW[23]=aryuWW[23]*aryuWW[21];
   aryuWW[23]=1./6.*aryuWW[23] - aryuWW[25] + 53./6.*aryuWW[15] + 3*
   aryuWW[14];
   aryuWW[23]=aryuWW[23]*aryuWW[21];
   aryuWW[23]=aryuWW[24] + aryuWW[23];
   aryuWW[18]=aryuWW[18] + 1./2.*aryuWW[23] + aryuWW[22];
   aryuWW[18]=aryuWW[4]*aryuWW[18];
   aryuWW[22]= - 209./8. - 4*aryuWW[10];
   aryuWW[23]=aryuWW[7]*aryuWW[10];
   aryuWW[24]= - aryuWW[7] + aryuWW[6];
   aryuWW[24]=aryuWW[1]*aryuWW[24];
   aryuWW[22]=aryuWW[24] + 4./3.*aryuWW[23] + aryuWW[11] + 1./9.*
   aryuWW[22] - 33./4.*aryuWW[12];
   aryuWW[19]=aryuWW[19]*aryuWW[22];
   aryuWW[22]=aryuWW[21] + 17;
   aryuWW[22]=aryuWW[12]*aryuWW[22];
   aryuWW[22]= - 1./2. + aryuWW[22];
   aryuWW[21]=aryuWW[22]*aryuWW[21];
   aryuWW[20]=aryuWW[17]*aryuWW[20];
   aryuWW[22]= - 1 + aryuWW[12];
   aryuWW[18]=3./4.*aryuWW[20] + aryuWW[18] + 1./12.*aryuWW[21] + 4*
   aryuWW[22] + aryuWW[19];

      yuWWret = aryuWW[18]*aryuWW[2];
      return yuWWret;
}
