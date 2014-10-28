#include <bb.hpp>
std::complex<long double>
bb::mygl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[27], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMH,-1);
    aryubbGL[4]=pow(MMW,-1);
    aryubbGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryubbGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryubbGL[7]=pow(MMt,-1);
    aryubbGL[8]=Tsil::I2(0,MMH,MMt,mu2);
    aryubbGL[9]=Tsil::I2(0,0,MMH,mu2);
    aryubbGL[10]=Tsil::I2(0,0,MMt,mu2);
    aryubbGL[11]=Tsil::B(MMH,MMH,MMH,mu2);
    aryubbGL[12]=Tsil::A(MMH,mu2);
    aryubbGL[13]=Tsil::A(MMt,mu2);
    aryubbGL[14]=Tsil::B(MMH,MMt,MMt,mu2);
    aryubbGL[15]=Tsil::B(MMt,MMt,MMH,mu2);
    aryubbGL[16]=std::real(Tsil::B(0,0,MMH,mu2));
    aryubbGL[17]=std::real(Tsil::B(0,0,MMt,mu2));
    aryubbGL[18]=Tsil::Aeps(MMH,mu2);
    aryubbGL[19]=Tsil::Aeps(MMt,mu2);
   aryubbGL[20]=9*aryubbGL[16];
   aryubbGL[21]=27*aryubbGL[11];
   aryubbGL[22]=aryubbGL[3]*aryubbGL[12];
   aryubbGL[22]=11*aryubbGL[22] + 15*aryubbGL[14] + 3./4.*aryubbGL[17]
    + aryubbGL[21] + 61./16. + aryubbGL[20];
   aryubbGL[23]=aryubbGL[13]*aryubbGL[3];
   aryubbGL[22]=1./2.*aryubbGL[22] + 9*aryubbGL[23];
   aryubbGL[22]=aryubbGL[13]*aryubbGL[22];
   aryubbGL[23]= - 31*aryubbGL[10] + 75*aryubbGL[18];
   aryubbGL[24]=9*aryubbGL[15];
   aryubbGL[25]=49./8. + aryubbGL[24];
   aryubbGL[25]=aryubbGL[12]*aryubbGL[25];
   aryubbGL[24]= - 7./4.*aryubbGL[14] - 85./8. + aryubbGL[24];
   aryubbGL[24]=MMH*aryubbGL[24];
   aryubbGL[22]=1./4.*aryubbGL[24] + aryubbGL[22] + 1./4.*aryubbGL[25]
    - 35./8.*aryubbGL[6] - 5./16.*aryubbGL[8] + 1./16.*aryubbGL[23] + 
   11*aryubbGL[19];
   aryubbGL[23]=1./2.*aryubbGL[10] - aryubbGL[18];
   aryubbGL[24]= - 5./2. - 3*aryubbGL[15];
   aryubbGL[24]=aryubbGL[12]*aryubbGL[24];
   aryubbGL[23]=3*aryubbGL[24] + 27./2.*aryubbGL[6] - 1./2.*aryubbGL[8]
    + 13*aryubbGL[23] - 33*aryubbGL[19];
   aryubbGL[23]=aryubbGL[3]*aryubbGL[23];
   aryubbGL[24]=463./8. + 7*aryubbGL[17];
   aryubbGL[23]=aryubbGL[23] + 19./4.*aryubbGL[14] + 1./16.*
   aryubbGL[24] - 9*aryubbGL[15];
   aryubbGL[24]=3*aryubbGL[15] - 13./8. - aryubbGL[17];
   aryubbGL[24]=1./2.*aryubbGL[24] - 2*aryubbGL[14];
   aryubbGL[24]=aryubbGL[3]*aryubbGL[24];
   aryubbGL[25]=pow(aryubbGL[3],2);
   aryubbGL[26]= - aryubbGL[13]*aryubbGL[25];
   aryubbGL[24]=aryubbGL[24] + 3./2.*aryubbGL[26];
   aryubbGL[24]=aryubbGL[13]*aryubbGL[24];
   aryubbGL[26]=3 - 1./2.*aryubbGL[17];
   aryubbGL[26]=1./2.*aryubbGL[26] - aryubbGL[14];
   aryubbGL[26]=aryubbGL[3]*aryubbGL[26];
   aryubbGL[25]= - aryubbGL[13]*aryubbGL[25]*aryubbGL[15];
   aryubbGL[25]=aryubbGL[26] + 6*aryubbGL[25];
   aryubbGL[25]=MMt*aryubbGL[25];
   aryubbGL[23]=3*aryubbGL[25] + 1./4.*aryubbGL[23] + 3*aryubbGL[24];
   aryubbGL[23]=MMt*aryubbGL[23];
   aryubbGL[22]=1./4.*aryubbGL[22] + aryubbGL[23];
   aryubbGL[22]=MMt*aryubbGL[22];
   aryubbGL[20]=aryubbGL[21] + 31 + aryubbGL[20];
   aryubbGL[20]=aryubbGL[12]*aryubbGL[20];
   aryubbGL[21]= - aryubbGL[12]*aryubbGL[7];
   aryubbGL[23]= - aryubbGL[13]*aryubbGL[7];
   aryubbGL[21]=aryubbGL[23] + aryubbGL[21] + 1./2. - 3*aryubbGL[14];
   aryubbGL[21]=aryubbGL[13]*aryubbGL[21];
   aryubbGL[23]= - 15*aryubbGL[5] - 17*aryubbGL[9];
   aryubbGL[20]=1./4.*aryubbGL[21] + 1./8.*aryubbGL[20] - 13./8.*
   aryubbGL[6] + 19./8.*aryubbGL[8] + 7./8.*aryubbGL[19] + 1./8.*
   aryubbGL[23] + 7*aryubbGL[18];
   aryubbGL[21]=1./2.*aryubbGL[6] - 3./2.*aryubbGL[8] + aryubbGL[9] + 1.
   /2.*aryubbGL[19];
   aryubbGL[21]=aryubbGL[7]*aryubbGL[21];
   aryubbGL[21]=1./8.*aryubbGL[21] + 27./32.*aryubbGL[11] - 1 + 9./32.*
   aryubbGL[16];
   aryubbGL[21]=MMH*aryubbGL[21];
   aryubbGL[20]=1./4.*aryubbGL[20] + aryubbGL[21];
   aryubbGL[20]=MMH*aryubbGL[20];
   aryubbGL[21]= - aryubbGL[12] + 7./8.*aryubbGL[13];
   aryubbGL[21]=aryubbGL[13]*aryubbGL[21];
   aryubbGL[23]=pow(aryubbGL[12],2);
   aryubbGL[21]= - 1./8.*aryubbGL[23] + aryubbGL[21];
   aryubbGL[20]=9./8.*aryubbGL[21] + aryubbGL[20];
   aryubbGL[20]=1./2.*aryubbGL[20] + aryubbGL[22];

      yubbGLret = aryubbGL[20]*pow(aryubbGL[4],2)*pow(aryubbGL[2],4)*
      aryubbGL[1];
      return yubbGLret;
}
