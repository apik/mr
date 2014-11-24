#include <HH.hpp>
std::complex<long double>
HH<OS>::y10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHH[28], yuHHret;

    aryuHH[1]=double(nH);
    aryuHH[2]=double(boson);
    aryuHH[3]=pow(CW,-1);
    aryuHH[4]=pow(MMZ,-1);
    aryuHH[5]=pow(SW,-1);
    aryuHH[6]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHH[7]=pow(MMH,-1);
    aryuHH[8]=Tsil::B(MMb,MMb,MMH,mu2);
    aryuHH[9]=Tsil::A(MMt,mu2);
    aryuHH[10]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHH[11]=Tsil::B(MMZ,MMZ,MMH,mu2);
    aryuHH[12]=Tsil::B(MMW,MMW,MMH,mu2);
    aryuHH[13]=Tsil::A(MMZ,mu2);
    aryuHH[14]=Tsil::A(MMW,mu2);
    aryuHH[15]=1/( - MMb + MMt);
    aryuHH[16]=Tsil::A(MMb,mu2);
    aryuHH[17]=1/( - MMW + MMH);
    aryuHH[18]=Tsil::A(MMH,mu2);
   aryuHH[19]=aryuHH[6] - 1./2.;
   aryuHH[19]=aryuHH[19]*MMt;
   aryuHH[20]=aryuHH[8]*MMb;
   aryuHH[19]=aryuHH[19] + aryuHH[20] - aryuHH[9];
   aryuHH[20]=3./2.*aryuHH[1];
   aryuHH[19]=aryuHH[19]*aryuHH[20];
   aryuHH[20]=1./4.*aryuHH[12] + 1./8.;
   aryuHH[20]=MMH*aryuHH[20];
   aryuHH[19]=aryuHH[19] + aryuHH[13] + 2*aryuHH[14] + aryuHH[20];
   aryuHH[20]=pow(aryuHH[5],2);
   aryuHH[21]= - aryuHH[13] + aryuHH[14];
   aryuHH[21]=aryuHH[21]*aryuHH[20];
   aryuHH[21]=3./4.*aryuHH[21] + aryuHH[19];
   aryuHH[21]=aryuHH[21]*aryuHH[20];
   aryuHH[22]=pow(aryuHH[3],2);
   aryuHH[19]=aryuHH[19]*aryuHH[22];
   aryuHH[23]=aryuHH[22] + aryuHH[20];
   aryuHH[24]=aryuHH[23]*MMH;
   aryuHH[25]=aryuHH[10]*aryuHH[24];
   aryuHH[26]=aryuHH[23]*aryuHH[1];
   aryuHH[27]=MMb - MMt;
   aryuHH[27]=aryuHH[9] - aryuHH[16] - 1./2.*aryuHH[27];
   aryuHH[27]= - aryuHH[15]*aryuHH[27]*MMb*aryuHH[26];
   aryuHH[19]=3./2.*aryuHH[27] + 9./8.*aryuHH[25] + aryuHH[21] + 
   aryuHH[19];
   aryuHH[19]=aryuHH[4]*aryuHH[19];
   aryuHH[21]=aryuHH[8]*pow(MMb,2);
   aryuHH[25]=aryuHH[6]*pow(MMt,2);
   aryuHH[21]=aryuHH[21] + aryuHH[25];
   aryuHH[21]= - aryuHH[4]*aryuHH[21]*aryuHH[26];
   aryuHH[25]= - 1 + aryuHH[12];
   aryuHH[25]=aryuHH[25]*aryuHH[20];
   aryuHH[25]=aryuHH[25] - aryuHH[12];
   aryuHH[26]=aryuHH[11]*aryuHH[23];
   aryuHH[25]=3./2.*aryuHH[26] - aryuHH[22] + 2 + 3*aryuHH[25];
   aryuHH[25]=MMZ*aryuHH[25];
   aryuHH[21]=6*aryuHH[21] + aryuHH[25];
   aryuHH[21]=aryuHH[7]*aryuHH[21];
   aryuHH[24]=aryuHH[4]*aryuHH[24];
   aryuHH[23]=1./4.*aryuHH[24] - aryuHH[23];
   aryuHH[23]=aryuHH[11]*aryuHH[23];
   aryuHH[24]=aryuHH[18] - aryuHH[14];
   aryuHH[24]=aryuHH[17]*aryuHH[24];
   aryuHH[24]=3./4.*aryuHH[24] + 3./8. - aryuHH[12];
   aryuHH[20]=aryuHH[24]*aryuHH[20];
   aryuHH[19]=aryuHH[21] + 1./2.*aryuHH[23] + aryuHH[20] + 1./8.*
   aryuHH[22] + aryuHH[19];

      yuHHret = aryuHH[19]*aryuHH[2];
      return yuHHret;
}
