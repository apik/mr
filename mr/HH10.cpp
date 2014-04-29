#include <HH.hpp>
std::complex<long double> HH::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHH[20];

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
   mHH[13]=pow(mHH[2],2);
   mHH[14]=pow(mHH[4],2);
   mHH[13]=mHH[13] + mHH[14];
   mHH[15]=mHH[13]*mHH[3];
   mHH[16]=1./4.*MMH;
   mHH[16]=mHH[15]*mHH[16];
   mHH[17]=3*mHH[6];
   mHH[18]=mHH[17]*MMZ;
   mHH[19]=mHH[18] - 1;
   mHH[19]=mHH[13]*mHH[19];
   mHH[19]=mHH[16] + mHH[19];
   mHH[19]=mHH[5]*mHH[19];
   mHH[20]=1./2.*mHH[15];
   mHH[13]=mHH[13]*mHH[17];
   mHH[13]=mHH[20] + mHH[13];
   mHH[13]=mHH[10]*mHH[13];
   mHH[13]=mHH[19] + mHH[13];
   mHH[17]=mHH[14]*mHH[17];
   mHH[17]=mHH[20] + mHH[17];
   mHH[17]=mHH[11]*mHH[17];
   mHH[19]=MMt*mHH[6];
   mHH[19]=2*mHH[19];
   mHH[19]=mHH[15]*mHH[19];
   mHH[20]=mHH[20] - mHH[19];
   mHH[20]=mHH[8]*MMt*mHH[20];
   mHH[19]= - mHH[12]*mHH[19];
   mHH[19]=mHH[19] + mHH[20];
   mHH[19]=nH*mHH[19];
   mHH[20]= - 1 + mHH[14];
   mHH[18]=mHH[20]*mHH[18];
   mHH[14]=mHH[18] - mHH[14] + mHH[16];
   mHH[14]=mHH[7]*mHH[14];
   mHH[16]=mHH[1]*MMH;
   mHH[16]=9./8.*mHH[16] + 3./4.*mHH[9];
   mHH[15]=mHH[15]*mHH[16];

      return 1./2.*mHH[13] + mHH[14] + mHH[15] + mHH[17] + 3*mHH[19];
}
