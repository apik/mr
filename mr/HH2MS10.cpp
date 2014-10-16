#include <HH.hpp>
std::complex<long double> hh::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHH2MS[21], mHH2MSret;

    mHH2MS[1]=double(nH);
    mHH2MS[2]=Tsil::B(mmt,mmt,mmH,mu2);
    mHH2MS[3]=pow(c,-1);
    mHH2MS[4]=pow(mmH,-1);
    mHH2MS[5]=pow(mmZ,-1);
    mHH2MS[6]=pow(s,-1);
    mHH2MS[7]=Tsil::A(mmt,mu2);
    mHH2MS[8]=Tsil::B(mmH,mmH,mmH,mu2);
    mHH2MS[9]=Tsil::B(mmZ,mmZ,mmH,mu2);
    mHH2MS[10]=Tsil::B(mmW,mmW,mmH,mu2);
    mHH2MS[11]=Tsil::A(mmH,mu2);
    mHH2MS[12]=Tsil::A(mmZ,mu2);
    mHH2MS[13]=Tsil::A(mmW,mu2);
   mHH2MS[14]=mmt*mHH2MS[2];
   mHH2MS[14]=mHH2MS[14] + mHH2MS[7];
   mHH2MS[15]=2*mHH2MS[4];
   mHH2MS[14]=mHH2MS[14]*mHH2MS[15];
   mHH2MS[14]=mHH2MS[14] - 1./2.*mHH2MS[2];
   mHH2MS[15]=mHH2MS[1]*mmt;
   mHH2MS[15]=3*mHH2MS[15];
   mHH2MS[14]=mHH2MS[14]*mHH2MS[15];
   mHH2MS[15]=mHH2MS[10] + 1./2.*mHH2MS[9];
   mHH2MS[16]=mHH2MS[15] + 9./2.*mHH2MS[8];
   mHH2MS[17]=1./4.*mmH;
   mHH2MS[16]=mHH2MS[16]*mHH2MS[17];
   mHH2MS[17]=mHH2MS[13] + 1./2.*mHH2MS[12];
   mHH2MS[18]=mHH2MS[17] + 3./2.*mHH2MS[11];
   mHH2MS[14]= - mHH2MS[14] + mHH2MS[16] + 1./2.*mHH2MS[18];
   mHH2MS[14]=mHH2MS[14]*mHH2MS[5];
   mHH2MS[16]=3*mHH2MS[4];
   mHH2MS[18]= - mHH2MS[12]*mHH2MS[16];
   mHH2MS[16]=mHH2MS[16]*mmZ;
   mHH2MS[19]=mHH2MS[16] - 1;
   mHH2MS[20]= - mHH2MS[9]*mHH2MS[19];
   mHH2MS[18]=mHH2MS[18] + mHH2MS[20];
   mHH2MS[18]=1./2.*mHH2MS[18] - mHH2MS[14];
   mHH2MS[18]=mHH2MS[18]*pow(mHH2MS[3],2);
   mHH2MS[15]= - mHH2MS[19]*mHH2MS[15];
   mHH2MS[17]= - mHH2MS[4]*mHH2MS[17];
   mHH2MS[14]= - mHH2MS[14] + 3*mHH2MS[17] + mHH2MS[15];
   mHH2MS[14]=mHH2MS[14]*pow(mHH2MS[6],2);
   mHH2MS[15]=mHH2MS[10]*mHH2MS[16];

      mHH2MSret = mHH2MS[14] + mHH2MS[15] + mHH2MS[18];
      return mHH2MSret;
}
