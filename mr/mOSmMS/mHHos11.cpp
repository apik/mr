#include <HH.hpp>
namespace mr
{
  long double HH<MS>::x11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armHHos[17], mHHosret;

    armHHos[1]=double(nH);
    armHHos[2]=pow(mmZ,-1);
    armHHos[3]=pow(mmH,-1);
    armHHos[4]=pow(s,-1);
    armHHos[5]=pow(c,-1);
    armHHos[6]=Tsil::B(mmt,mmt,mmH,mu2);
    armHHos[7]=Tsil::A(mmt,mu2);
    armHHos[8]=Tsil::Beps(mmt,mmt,mmH,mu2);
    armHHos[9]=Tsil::Aeps(mmt,mu2);
    armHHos[10]=prottttt0->M(0);
    armHHos[11]=prottttt0->Vxzuv(0);
    armHHos[12]=prottttt0->Suxv(0);
    armHHos[13]=2*armHHos[11];
    armHHos[14]=armHHos[13] + armHHos[10];
    armHHos[15]=8*mmt;
    armHHos[14]=armHHos[14]*armHHos[15];
    armHHos[15]=20 + 7*armHHos[6];
    armHHos[15]=armHHos[15]*armHHos[6];
    armHHos[14]= - armHHos[15] + armHHos[14] - 11 + 4*armHHos[8];
    armHHos[14]=armHHos[14]*mmt;
    armHHos[15]= - 1 + 3*armHHos[6];
    armHHos[16]=2*armHHos[7];
    armHHos[15]=armHHos[15]*armHHos[16];
    armHHos[15]=armHHos[15] + 3*armHHos[9];
    armHHos[14]= - armHHos[12] + armHHos[14] + 2*armHHos[15];
    armHHos[14]=armHHos[14]*mmt;
    armHHos[15]=pow(armHHos[7],2);
    armHHos[14]=armHHos[14] + 18*armHHos[15];
    armHHos[15]=armHHos[3]*armHHos[1];
    armHHos[15]=2*armHHos[15];
    armHHos[14]=armHHos[14]*armHHos[15];
    armHHos[13]=armHHos[13] + 3*armHHos[10];
    armHHos[15]=4*mmt;
    armHHos[13]=armHHos[13]*armHHos[15];
    armHHos[15]=armHHos[6] + 2;
    armHHos[15]=armHHos[15]*armHHos[6];
    armHHos[15]=armHHos[15] + 4;
    armHHos[16]=armHHos[10]*mmH;
    armHHos[13]=2*armHHos[16] - armHHos[13] + 3*armHHos[15];
    armHHos[13]=armHHos[1]*mmt*armHHos[13];
    armHHos[13]=armHHos[14] + armHHos[13];
    armHHos[14]= - pow(armHHos[4],2);
    armHHos[15]= - pow(armHHos[5],2);
    armHHos[14]=armHHos[14] + armHHos[15];

    mHHosret = 2*armHHos[14]*armHHos[13]*armHHos[2];
    return mHHosret.real();
  }
} // namespace mr
