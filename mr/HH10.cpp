#include <HH.hpp>
std::complex<long double> HH::m10(size_t nG)
{     
      
      
    std::complex<long double> mHH[22];

    mHH[1]=Tsil::B(MMH,MMH,MMH,mu2);
    mHH[2]=pow(CW,-1);
    mHH[3]=pow(MMZ,-1);
    mHH[4]=pow(SW,-1);
    mHH[5]=Tsil::B(MMZ,MMZ,MMH,mu2);
    mHH[6]=pow(MMH,-1);
    mHH[7]=Tsil::B(MMW,MMW,MMH,mu2);
    mHH[8]=Tsil::B(MMt,MMt,MMH,mu2);
    mHH[9]=Tsil::A(MMH,mu2);
    mHH[10]=Tsil::A(MMZ,mu2);
    mHH[11]=Tsil::A(MMW,mu2);
    mHH[12]=Tsil::A(MMt,mu2);
   mHH[13]=1./2.*mHH[5];
   mHH[14]=mHH[13] + mHH[7];
   mHH[15]=MMH*mHH[14];
   mHH[15]=3*mHH[9] + mHH[15];
   mHH[16]=mHH[11] + 1./2.*mHH[10];
   mHH[17]=mHH[1]*MMH;
   mHH[18]=MMt*mHH[8];
   mHH[15]=mHH[16] + 9./4.*mHH[17] + 3*mHH[18] + 1./2.*mHH[15];
   mHH[17]=pow(mHH[4],2);
   mHH[19]=pow(mHH[2],2);
   mHH[20]=mHH[17] + mHH[19];
   mHH[20]=mHH[20]*mHH[3];
   mHH[15]=mHH[15]*mHH[20];
   mHH[21]=mHH[7]*MMZ;
   mHH[22]=MMZ*mHH[13];
   mHH[16]=mHH[21] + mHH[22] + mHH[16];
   mHH[16]=mHH[16]*mHH[17];
   mHH[22]=mHH[5]*MMZ;
   mHH[22]=mHH[22] + mHH[10];
   mHH[22]=mHH[22]*mHH[19];
   mHH[18]=mHH[18] + mHH[12];
   mHH[18]= - MMt*mHH[18]*mHH[20];
   mHH[16]=2*mHH[18] + 1./2.*mHH[22] - mHH[21] + mHH[16];
   mHH[16]=mHH[6]*mHH[16];
   mHH[14]= - mHH[14]*mHH[17];
   mHH[13]= - mHH[19]*mHH[13];

      return mHH[13] + mHH[14] + 1./2.*mHH[15] + 3*mHH[16];
}
