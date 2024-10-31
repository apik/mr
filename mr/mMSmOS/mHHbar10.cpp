#include <HH.hpp>
namespace mr
{
  double HH<OS>::x10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> armHHbar[26], mHHbarret;

    armHHbar[1]=double(nH);
    armHHbar[2]=double(boson);
    armHHbar[3]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbar[4]=pow(CW,-1);
    armHHbar[5]=pow(MMH,-1);
    armHHbar[6]=pow(MMZ,-1);
    armHHbar[7]=pow(SW,-1);
    armHHbar[8]=Tsil::B(MMb,MMb,MMH,mu2);
    armHHbar[9]=Tsil::A(MMt,mu2);
    armHHbar[10]=Tsil::A(MMb,mu2);
    armHHbar[11]=Tsil::B(MMH,MMH,MMH,mu2);
    armHHbar[12]=Tsil::B(MMZ,MMZ,MMH,mu2);
    armHHbar[13]=Tsil::B(MMW,MMW,MMH,mu2);
    armHHbar[14]=Tsil::A(MMH,mu2);
    armHHbar[15]=Tsil::A(MMZ,mu2);
    armHHbar[16]=Tsil::A(MMW,mu2);
    armHHbar[17]=pow(armHHbar[4],2);
    armHHbar[18]=pow(armHHbar[7],2);
    armHHbar[17]=armHHbar[17] + armHHbar[18];
    armHHbar[19]=MMb*armHHbar[17]*armHHbar[8];
    armHHbar[20]= - armHHbar[10]*armHHbar[17];
    armHHbar[20]= - armHHbar[19] + armHHbar[20];
    armHHbar[20]=MMb*armHHbar[20];
    armHHbar[21]=MMt*armHHbar[3];
    armHHbar[22]= - armHHbar[21] - armHHbar[9];
    armHHbar[22]=MMt*armHHbar[17]*armHHbar[22];
    armHHbar[20]=armHHbar[22] + armHHbar[20];
    armHHbar[20]=armHHbar[6]*armHHbar[1]*armHHbar[20];
    armHHbar[22]=1./2.*armHHbar[15];
    armHHbar[22]=armHHbar[17]*armHHbar[22];
    armHHbar[23]=1./2.*armHHbar[12];
    armHHbar[23]=armHHbar[17]*armHHbar[23];
    armHHbar[24]= - 1 + armHHbar[18];
    armHHbar[24]=armHHbar[13]*armHHbar[24];
    armHHbar[24]=armHHbar[24] + armHHbar[23];
    armHHbar[24]=MMZ*armHHbar[24];
    armHHbar[25]=armHHbar[16]*armHHbar[18];
    armHHbar[20]=armHHbar[24] + 2*armHHbar[20] + armHHbar[25] + 
      armHHbar[22];
    armHHbar[24]=3*armHHbar[5];
    armHHbar[20]=armHHbar[24]*armHHbar[20];
    armHHbar[21]=armHHbar[17]*armHHbar[21];
    armHHbar[19]=armHHbar[21] + armHHbar[19];
    armHHbar[19]=armHHbar[1]*armHHbar[19];
    armHHbar[21]=9./4.*armHHbar[11] + 1./4.*armHHbar[12] + 1./2.*
      armHHbar[13];
    armHHbar[21]=armHHbar[21]*MMH;
    armHHbar[21]=armHHbar[21] + armHHbar[16] + 3./2.*armHHbar[14];
    armHHbar[17]=armHHbar[17]*armHHbar[21];
    armHHbar[17]=armHHbar[22] + 3*armHHbar[19] + armHHbar[17];
    armHHbar[17]=armHHbar[6]*armHHbar[17];
    armHHbar[18]= - armHHbar[13]*armHHbar[18];
    armHHbar[17]=1./2.*armHHbar[17] + armHHbar[18] - armHHbar[23] + 
      armHHbar[20];

    mHHbarret = armHHbar[17]*armHHbar[2];
    return mHHbarret.real();
  }
} // namespace mr
