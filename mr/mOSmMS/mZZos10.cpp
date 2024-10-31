#include <ZZ.hpp>
namespace mr
{
  double ZZ<MS>::x10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> armZZos[27], mZZosret;

    armZZos[1]=double(nL + nH);
    armZZos[2]=pow(s,-1);
    armZZos[3]=pow(c,-1);
    armZZos[4]=std::real(Tsil::B(0,0,mmZ,mu2));
    armZZos[5]=double(nH);
    armZZos[6]=pow(mmZ,-1);
    armZZos[7]=Tsil::B(mmt,mmt,mmZ,mu2);
    armZZos[8]=Tsil::A(mmt,mu2);
    armZZos[9]=pow(mmH,-1);
    armZZos[10]=double(nL);
    armZZos[11]=double(boson);
    armZZos[12]=Tsil::B(mmW,mmW,mmZ,mu2);
    armZZos[13]=Tsil::B(mmZ,mmH,mmZ,mu2);
    armZZos[14]=Tsil::A(mmW,mu2);
    armZZos[15]=Tsil::A(mmZ,mu2);
    armZZos[16]=Tsil::A(mmH,mu2);
    armZZos[17]=pow(armZZos[3],2);
    armZZos[18]=pow(armZZos[2],2);
    armZZos[19]=armZZos[17] + armZZos[18];
    armZZos[20]= - armZZos[15]*armZZos[19];
    armZZos[21]=3*armZZos[18];
    armZZos[22]= - armZZos[17] + 2 - armZZos[21];
    armZZos[22]=mmZ*armZZos[22];
    armZZos[20]=armZZos[22] + 3./2.*armZZos[20];
    armZZos[20]=armZZos[9]*armZZos[20];
    armZZos[22]=armZZos[19]*armZZos[6];
    armZZos[23]=armZZos[15] - armZZos[16];
    armZZos[23]=armZZos[23]*armZZos[22];
    armZZos[23]=1./6.*armZZos[23] + 1./3.*armZZos[19];
    armZZos[23]=mmH*armZZos[23];
    armZZos[24]=armZZos[16] + 1./3.*armZZos[15];
    armZZos[24]= - armZZos[24]*armZZos[19];
    armZZos[23]=armZZos[23] + armZZos[24];
    armZZos[23]=armZZos[6]*armZZos[23];
    armZZos[24]=armZZos[6]*pow(mmH,2);
    armZZos[24]=mmH - 1./4.*armZZos[24];
    armZZos[22]=armZZos[24]*armZZos[22];
    armZZos[22]=1./3.*armZZos[22] - armZZos[19];
    armZZos[22]=armZZos[13]*armZZos[22];
    armZZos[24]=1./6.*armZZos[17];
    armZZos[25]=armZZos[24] - 8 + 59./6.*armZZos[18];
    armZZos[26]= - 1./12.*armZZos[17] - 29./3. + 33./4.*armZZos[18];
    armZZos[26]=armZZos[12]*armZZos[26];
    armZZos[24]= - armZZos[24] - 4 + 5./2.*armZZos[18];
    armZZos[24]=armZZos[6]*armZZos[24];
    armZZos[21]= - armZZos[9]*armZZos[21];
    armZZos[21]=armZZos[21] + armZZos[24];
    armZZos[21]=armZZos[14]*armZZos[21];
    armZZos[24]= - 1 - armZZos[12];
    armZZos[24]=armZZos[24]*pow(c,2);
    armZZos[20]=armZZos[21] + 4*armZZos[24] + armZZos[26] + armZZos[22]
      + 1./2.*armZZos[23] + 1./3.*armZZos[25] + armZZos[20];
    armZZos[20]=armZZos[11]*armZZos[20];
    armZZos[21]=17./9.*armZZos[17] + armZZos[18] - 32./9.;
    armZZos[19]=armZZos[8]*armZZos[9]*armZZos[19];
    armZZos[22]=1./2.*armZZos[18];
    armZZos[23]= - 7./18.*armZZos[17] + 32./9. + armZZos[22];
    armZZos[23]=armZZos[7]*armZZos[23];
    armZZos[19]=armZZos[23] + 6*armZZos[19] - armZZos[21];
    armZZos[19]=mmt*armZZos[19];
    armZZos[21]= - armZZos[8]*armZZos[21];
    armZZos[19]=armZZos[19] + armZZos[21];
    armZZos[19]=armZZos[19]*armZZos[6];
    armZZos[21]=11./9.*armZZos[17] + armZZos[18] - 20./9.;
    armZZos[23]= - 5./18.*armZZos[17] + 4./9. - armZZos[22];
    armZZos[23]=armZZos[4]*armZZos[23];
    armZZos[22]= - 17./18.*armZZos[17] + 16./9. - armZZos[22];
    armZZos[22]=armZZos[7]*armZZos[22];
    armZZos[19]=armZZos[19] + armZZos[22] + 1./3.*armZZos[21] + 
      armZZos[23];
    armZZos[19]=armZZos[5]*armZZos[19];
    armZZos[21]= - armZZos[10]*armZZos[21];
    armZZos[18]=armZZos[18] - 4;
    armZZos[17]=armZZos[17] + 1./3.*armZZos[18];
    armZZos[17]= - armZZos[1]*armZZos[17];
    armZZos[17]=armZZos[17] + armZZos[21];
    armZZos[18]=armZZos[4] - 1./3.;
    armZZos[17]=armZZos[18]*armZZos[17];

    mZZosret = armZZos[17] + armZZos[19] + armZZos[20];
    return mZZosret.real();
  }
} // namespace mr
