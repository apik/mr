#include <HH.hpp>
std::complex<long double>
HH<MS>::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHos[19], mHHosret;

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
   armHHos[13]=2*armHHos[11] + 3*armHHos[10];
   armHHos[14]=pow(armHHos[5],2);
   armHHos[15]=armHHos[14]*armHHos[13];
   armHHos[16]=pow(armHHos[4],2);
   armHHos[13]=armHHos[16]*armHHos[13];
   armHHos[13]=armHHos[15] + armHHos[13];
   armHHos[15]=20 + 7*armHHos[6];
   armHHos[15]=armHHos[6]*armHHos[15];
   armHHos[15]=armHHos[15] + 11 - 4*armHHos[8];
   armHHos[17]=armHHos[14]*armHHos[15];
   armHHos[15]=armHHos[16]*armHHos[15];
   armHHos[15]=armHHos[17] + armHHos[15];
   armHHos[15]=armHHos[3]*armHHos[15];
   armHHos[17]= - 2*armHHos[11] - armHHos[10];
   armHHos[18]=armHHos[14]*armHHos[17];
   armHHos[17]=armHHos[16]*armHHos[17];
   armHHos[17]=armHHos[18] + armHHos[17];
   armHHos[17]=mmt*armHHos[3]*armHHos[17];
   armHHos[13]=8*armHHos[17] + 2*armHHos[13] + armHHos[15];
   armHHos[13]=mmt*armHHos[13];
   armHHos[15]= - armHHos[6]*armHHos[7];
   armHHos[15]=12*armHHos[15] + 4*armHHos[7] + armHHos[12] - 6*
   armHHos[9];
   armHHos[17]=armHHos[14]*armHHos[15];
   armHHos[15]=armHHos[16]*armHHos[15];
   armHHos[15]=armHHos[17] + armHHos[15];
   armHHos[15]=armHHos[3]*armHHos[15];
   armHHos[17]= - armHHos[10]*mmH;
   armHHos[17]= - 6 + armHHos[17];
   armHHos[18]= - 2 - armHHos[6];
   armHHos[18]=armHHos[6]*armHHos[18];
   armHHos[17]=2*armHHos[17] + 3*armHHos[18];
   armHHos[18]=armHHos[14]*armHHos[17];
   armHHos[17]=armHHos[16]*armHHos[17];
   armHHos[13]=2*armHHos[13] + 2*armHHos[15] + armHHos[18] + 
   armHHos[17];
   armHHos[13]=mmt*armHHos[13];
   armHHos[15]=pow(armHHos[7],2);
   armHHos[14]= - armHHos[14]*armHHos[15];
   armHHos[15]= - armHHos[16]*armHHos[15];
   armHHos[14]=armHHos[14] + armHHos[15];
   armHHos[14]=armHHos[3]*armHHos[14];
   armHHos[13]=36*armHHos[14] + armHHos[13];

      mHHosret = 2*armHHos[13]*armHHos[2]*armHHos[1];
      return mHHosret;
}
