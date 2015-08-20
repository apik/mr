#include <WW.hpp>
namespace mr
{
  long double WW<OS>::x10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armWWbar[26], mWWbarret;

    armWWbar[1]=double(nH);
    armWWbar[2]=double(boson);
    armWWbar[3]=pow(CW,-1);
    armWWbar[4]=pow(MMZ,-1);
    armWWbar[5]=pow(SW,-1);
    armWWbar[6]=Tsil::B(MMt,MMb,MMW,mu2);
    armWWbar[7]=Tsil::B(0,0,MMW,mu2);
    armWWbar[8]=Tsil::A(MMt,mu2);
    armWWbar[9]=pow(MMH,-1);
    armWWbar[10]=Tsil::A(MMb,mu2);
    armWWbar[11]=double(nL + nH);
    armWWbar[12]=Tsil::B(MMW,MMH,MMW,mu2);
    armWWbar[13]=Tsil::B(MMW,MMZ,MMW,mu2);
    armWWbar[14]=Tsil::A(MMH,mu2);
    armWWbar[15]=Tsil::A(MMZ,mu2);
    armWWbar[16]=Tsil::A(MMW,mu2);
    armWWbar[17]=6*armWWbar[9];
    armWWbar[18]= - armWWbar[10]*armWWbar[17];
    armWWbar[18]=armWWbar[18] + 1;
    armWWbar[19]=MMb*armWWbar[1];
    armWWbar[18]=armWWbar[19]*armWWbar[18];
    armWWbar[20]=armWWbar[1]*MMt;
    armWWbar[21]= - armWWbar[20] - armWWbar[19];
    armWWbar[21]=armWWbar[6]*armWWbar[21];
    armWWbar[17]= - armWWbar[20]*armWWbar[17];
    armWWbar[17]=armWWbar[1] + armWWbar[17];
    armWWbar[17]=armWWbar[8]*armWWbar[17];
    armWWbar[22]=armWWbar[10]*armWWbar[1];
    armWWbar[23]= - 1./2. - armWWbar[12];
    armWWbar[23]=MMH*armWWbar[23];
    armWWbar[17]=armWWbar[17] + 1./3.*armWWbar[23] + 1./2.*armWWbar[21]
      + armWWbar[18] + armWWbar[22] + 1./2.*armWWbar[14] + armWWbar[20];
    armWWbar[18]=pow(armWWbar[5],2);
    armWWbar[21]=pow(armWWbar[3],2);
    armWWbar[22]=armWWbar[18] + armWWbar[21];
    armWWbar[17]=armWWbar[22]*armWWbar[17];
    armWWbar[23]=armWWbar[21] + 1;
    armWWbar[23]=armWWbar[23]*armWWbar[21];
    armWWbar[23]=armWWbar[23] + armWWbar[18];
    armWWbar[19]=armWWbar[23]*armWWbar[19];
    armWWbar[20]=armWWbar[23]*armWWbar[20];
    armWWbar[24]=armWWbar[19] - armWWbar[20];
    armWWbar[25]=armWWbar[8] - armWWbar[10];
    armWWbar[25]=1./2.*armWWbar[25];
    armWWbar[24]=armWWbar[24]*armWWbar[25];
    armWWbar[19]=armWWbar[20] - 1./2.*armWWbar[19];
    armWWbar[19]=MMb*armWWbar[19];
    armWWbar[20]=1./2.*armWWbar[1];
    armWWbar[20]= - armWWbar[23]*pow(MMt,2)*armWWbar[20];
    armWWbar[19]=armWWbar[20] + armWWbar[19];
    armWWbar[19]=armWWbar[6]*armWWbar[19];
    armWWbar[20]=armWWbar[14]*armWWbar[23];
    armWWbar[23]=armWWbar[23]*MMH;
    armWWbar[25]=armWWbar[12]*armWWbar[23];
    armWWbar[20]=armWWbar[20] + armWWbar[25];
    armWWbar[20]=MMH*armWWbar[20];
    armWWbar[19]=1./12.*armWWbar[20] + armWWbar[19] + armWWbar[24];
    armWWbar[19]=armWWbar[4]*armWWbar[19];
    armWWbar[20]=3 + 1./3.*armWWbar[21];
    armWWbar[20]=armWWbar[20]*armWWbar[21];
    armWWbar[20]= - 5*armWWbar[18] + armWWbar[20];
    armWWbar[20]=armWWbar[15]*armWWbar[20];
    armWWbar[17]=armWWbar[19] + 1./4.*armWWbar[20] + armWWbar[17];
    armWWbar[17]=armWWbar[4]*armWWbar[17];
    armWWbar[19]=armWWbar[15]*armWWbar[22];
    armWWbar[20]=3*armWWbar[18];
    armWWbar[22]=armWWbar[21] - 2 + armWWbar[20];
    armWWbar[22]=MMZ*armWWbar[22];
    armWWbar[19]=3./2.*armWWbar[19] + armWWbar[22];
    armWWbar[19]=armWWbar[9]*armWWbar[19];
    armWWbar[22]=35 - armWWbar[21];
    armWWbar[22]=armWWbar[22]*armWWbar[21];
    armWWbar[23]= - armWWbar[4]*armWWbar[23];
    armWWbar[22]=armWWbar[23] - 13*armWWbar[18] + armWWbar[22];
    armWWbar[22]=armWWbar[4]*armWWbar[22];
    armWWbar[20]=armWWbar[9]*armWWbar[20];
    armWWbar[20]=armWWbar[20] + 1./12.*armWWbar[22];
    armWWbar[20]=armWWbar[16]*armWWbar[20];
    armWWbar[22]=armWWbar[6]*armWWbar[1];
    armWWbar[23]= - armWWbar[1] + 4./3.*armWWbar[11];
    armWWbar[23]=armWWbar[7]*armWWbar[23];
    armWWbar[22]=armWWbar[23] - 4./9.*armWWbar[11] + armWWbar[22] + 
      armWWbar[12] - 59./18.;
    armWWbar[22]=armWWbar[18]*armWWbar[22];
    armWWbar[23]=17 + armWWbar[21];
    armWWbar[23]=armWWbar[23]*armWWbar[21];
    armWWbar[18]=1./12.*armWWbar[23] + 4 - 33./4.*armWWbar[18];
    armWWbar[18]=armWWbar[13]*armWWbar[18];
    armWWbar[17]=armWWbar[20] + armWWbar[17] + armWWbar[18] + 
      armWWbar[19] - 1./6.*armWWbar[21] - 4 + armWWbar[22];

    mWWbarret = armWWbar[17]*armWWbar[2];
    return mWWbarret.real();
  }
} // namespace mr
