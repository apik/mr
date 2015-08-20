#include <HH.hpp>
namespace mr
{
  long double HH<MS>::x10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armHHos[22], mHHosret;

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
    armHHos[15]=2*armHHos[4];
    armHHos[16]=armHHos[2]*armHHos[15]*mmt;
    armHHos[15]=armHHos[15]*armHHos[7];
    armHHos[15]= - armHHos[16] - armHHos[15] + 1./2.*armHHos[2];
    armHHos[16]=pow(armHHos[5],2);
    armHHos[17]=pow(armHHos[6],2);
    armHHos[18]= - armHHos[16] - armHHos[17];
    armHHos[15]=armHHos[18]*armHHos[1]*mmt*armHHos[3]*armHHos[15];
    armHHos[18]=armHHos[9] + 9./2.*armHHos[11];
    armHHos[18]=1./4.*armHHos[10] + 1./2.*armHHos[18];
    armHHos[18]=armHHos[18]*mmH;
    armHHos[19]=armHHos[13] + 3*armHHos[14];
    armHHos[18]=armHHos[18] + armHHos[12] + 1./2.*armHHos[19];
    armHHos[18]=armHHos[18]*armHHos[3];
    armHHos[19]=3*armHHos[4];
    armHHos[20]=armHHos[19]*mmZ;
    armHHos[21]=armHHos[20] - 1;
    armHHos[21]=armHHos[21]*armHHos[10];
    armHHos[18]=armHHos[18] + armHHos[21];
    armHHos[21]=mmZ*armHHos[9];
    armHHos[21]= - armHHos[12] - 1./2.*armHHos[13] - armHHos[21];
    armHHos[21]=armHHos[21]*armHHos[19];
    armHHos[21]=armHHos[9] + armHHos[21] - 1./2.*armHHos[18];
    armHHos[16]=armHHos[21]*armHHos[16];
    armHHos[19]= - armHHos[13]*armHHos[19];
    armHHos[18]=armHHos[19] - armHHos[18];
    armHHos[17]=armHHos[18]*armHHos[17];
    armHHos[18]=armHHos[9]*armHHos[20];
    armHHos[16]=1./2.*armHHos[17] + armHHos[18] + armHHos[16];
    armHHos[16]=armHHos[8]*armHHos[16];

    mHHosret = 3*armHHos[15] + armHHos[16];
    return mHHosret.real();
  }
} // namespace mr
