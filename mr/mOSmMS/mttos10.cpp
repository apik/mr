#include <tt.hpp>
namespace mr
{
  long double tt<MS>::x10(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armttos[24], mttosret;

    armttos[1]=double(nH);
    armttos[2]=Tsil::A(mmt,mu2);
    armttos[3]=pow(mmZ,-1);
    armttos[4]=pow(mmH,-1);
    armttos[5]=pow(s,-1);
    armttos[6]=pow(c,-1);
    armttos[7]=double(boson);
    armttos[8]=Tsil::B(mmZ,mmt,mmt,mu2);
    armttos[9]=pow(mmt,-1);
    armttos[10]=Tsil::B(mmH,mmt,mmt,mu2);
    armttos[11]=Tsil::A(mmW,mu2);
    armttos[12]=Tsil::A(mmZ,mu2);
    armttos[13]=Tsil::A(mmH,mu2);
    armttos[14]=std::real(Tsil::B(0,mmW,mmt,mu2));
    armttos[15]=1./4.*armttos[14];
    armttos[16]=pow(armttos[5],2);
    armttos[17]= - 1 + armttos[16];
    armttos[17]=armttos[17]*armttos[15];
    armttos[18]=1./8.*armttos[16];
    armttos[19]=pow(armttos[6],2);
    armttos[20]=armttos[18] + 17./72.*armttos[19];
    armttos[21]=armttos[20] - 4./9.;
    armttos[22]=armttos[8]*armttos[21];
    armttos[17]=armttos[17] + armttos[22];
    armttos[17]=mmZ*armttos[17];
    armttos[21]=armttos[12]*armttos[21];
    armttos[20]= - 8./9. - armttos[20];
    armttos[20]=armttos[2]*armttos[20];
    armttos[22]=1./4.*armttos[11];
    armttos[23]=armttos[16]*armttos[22];
    armttos[17]=armttos[23] + armttos[17] + armttos[21] + armttos[20];
    armttos[17]=armttos[9]*armttos[17];
    armttos[20]=armttos[13] + armttos[2];
    armttos[21]=armttos[10]*mmt;
    armttos[23]=mmH*armttos[10];
    armttos[20]=1./4.*armttos[23] - armttos[22] - armttos[21] - 1./2.*
      armttos[20];
    armttos[21]=armttos[19] + armttos[16];
    armttos[20]=armttos[21]*armttos[20];
    armttos[22]=armttos[21]*mmt;
    armttos[15]= - armttos[15]*armttos[22];
    armttos[15]=armttos[15] + armttos[20];
    armttos[15]=armttos[3]*armttos[15];
    armttos[20]=3./2.*armttos[4];
    armttos[23]= - armttos[11]*armttos[20];
    armttos[23]=armttos[23] - 1./8.*armttos[14] + 3./8.;
    armttos[23]=armttos[16]*armttos[23];
    armttos[21]= - armttos[12]*armttos[21]*armttos[4];
    armttos[18]=7./72.*armttos[19] - armttos[18] - 8./9.;
    armttos[18]=armttos[8]*armttos[18];
    armttos[16]= - armttos[16]*armttos[20];
    armttos[20]= - 1./2.*armttos[19] + 1;
    armttos[20]=armttos[4]*armttos[20];
    armttos[16]=armttos[16] + armttos[20];
    armttos[16]=mmZ*armttos[16];
    armttos[15]=1./2.*armttos[15] + armttos[17] + armttos[16] + 
      armttos[18] + 3./4.*armttos[21] + 1./72.*armttos[19] + 8./9. + 
      armttos[23];
    armttos[15]=armttos[7]*armttos[15];
    armttos[16]=armttos[3]*armttos[2]*armttos[1]*armttos[4]*armttos[22];

    mttosret = armttos[15] + 3*armttos[16];
    return mttosret.real();
  }
} // namespace mr
