#include <HH.hpp>
std::complex<long double>
HH<MS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHos[21], mHHosret;

    armHHos[1]=double(nH);
    armHHos[2]=Tsil::B(mmt,mmt,mmH,mu2);
    armHHos[3]=pow(mmZ,-1);
    armHHos[4]=pow(mmH,-1);
    armHHos[5]=pow(s,-1);
    armHHos[6]=pow(c,-1);
    armHHos[7]=Tsil::A(mmt,mu2);
    armHHos[8]=double(boson);
    armHHos[9]=Tsil::B(mmW,mmW,mmH,mu2);
    armHHos[10]=Tsil::B(mmZ,mmZ,mmH,mu2);
    armHHos[11]=Tsil::B(mmH,mmH,mmH,mu2);
    armHHos[12]=Tsil::A(mmW,mu2);
    armHHos[13]=Tsil::A(mmZ,mu2);
    armHHos[14]=Tsil::A(mmH,mu2);
   armHHos[15]= - 1./2.*armHHos[13];
   armHHos[16]= - 9./2.*armHHos[11] - armHHos[9];
   armHHos[16]=mmH*armHHos[16];
   armHHos[17]= - armHHos[10]*mmH;
   armHHos[16]=1./4.*armHHos[17] + 1./2.*armHHos[16] + armHHos[15] - 3./
   2.*armHHos[14] - armHHos[12];
   armHHos[17]=pow(armHHos[6],2);
   armHHos[18]=armHHos[17]*armHHos[16];
   armHHos[19]=pow(armHHos[5],2);
   armHHos[16]=armHHos[19]*armHHos[16];
   armHHos[16]=armHHos[18] + armHHos[16];
   armHHos[16]=armHHos[3]*armHHos[16];
   armHHos[18]= - armHHos[9]*mmZ;
   armHHos[20]= - armHHos[10]*mmZ;
   armHHos[15]=1./2.*armHHos[20] + armHHos[18] - armHHos[12] + 
   armHHos[15];
   armHHos[15]=armHHos[4]*armHHos[15];
   armHHos[15]=3*armHHos[15] + armHHos[9] + 1./2.*armHHos[10];
   armHHos[15]=armHHos[19]*armHHos[15];
   armHHos[18]= - armHHos[13] + armHHos[20];
   armHHos[18]=armHHos[4]*armHHos[18];
   armHHos[18]=armHHos[10] + 3*armHHos[18];
   armHHos[18]=armHHos[17]*armHHos[18];
   armHHos[20]=armHHos[4]*armHHos[9]*mmZ;
   armHHos[15]=1./2.*armHHos[16] + armHHos[15] + 3*armHHos[20] + 1./2.*
   armHHos[18];
   armHHos[15]=armHHos[8]*armHHos[15];
   armHHos[16]= - armHHos[1]*mmt*armHHos[2];
   armHHos[18]=mmt*armHHos[2];
   armHHos[18]=armHHos[7] + armHHos[18];
   armHHos[18]=armHHos[4]*armHHos[1]*mmt*armHHos[18];
   armHHos[16]=1./2.*armHHos[16] + 2*armHHos[18];
   armHHos[17]=armHHos[17]*armHHos[16];
   armHHos[16]=armHHos[19]*armHHos[16];
   armHHos[16]=armHHos[17] + armHHos[16];
   armHHos[16]=armHHos[3]*armHHos[16];

      mHHosret = armHHos[15] + 3*armHHos[16];
      return mHHosret;
}
