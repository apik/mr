#include <HH.hpp>
std::complex<long double> HH::m10(size_t nG)
{     
      
      
    std::complex<long double> mHH[18];

    mHH[1]=Tsil::B(MMH,MMH,MMH,mu2);
    mHH[2]=pow(MMW,-1);
    mHH[3]=Tsil::B(MMZ,MMZ,MMH,mu2);
    mHH[4]=pow(MMH,-1);
    mHH[5]=Tsil::B(MMW,MMW,MMH,mu2);
    mHH[6]=Tsil::B(MMt,MMt,MMH,mu2);
    mHH[7]=Tsil::A(MMH,mu2);
    mHH[8]=Tsil::A(MMZ,mu2);
    mHH[9]=Tsil::A(MMW,mu2);
    mHH[10]=Tsil::A(MMt,mu2);
    mHH[11]=1/(1 - pow(MMZ,-1)*MMW);
    mHH[12]=pow(MMZ,-1);
   mHH[13]=9./2.*mHH[1] + mHH[5];
   mHH[13]=MMH*mHH[13];
   mHH[14]=MMt*mHH[6];
   mHH[15]=1./2.*mHH[8];
   mHH[13]=3*mHH[14] + 1./2.*mHH[13] + mHH[15] + 3./2.*mHH[7] + mHH[9];
   mHH[14]=1./4.*MMH - MMZ;
   mHH[14]=mHH[3]*mHH[14];
   mHH[14]=mHH[13] + mHH[14];
   mHH[16]=MMZ*mHH[8];
   mHH[17]=mHH[3]*pow(MMZ,2);
   mHH[18]= - MMt*mHH[6];
   mHH[18]= - mHH[10] + mHH[18];
   mHH[18]=MMt*mHH[18];
   mHH[16]=1./2.*mHH[17] + 1./2.*mHH[16] + 2*mHH[18];
   mHH[16]=mHH[4]*mHH[16];
   mHH[14]=1./2.*mHH[14] + 3*mHH[16];
   mHH[14]=mHH[2]*mHH[14];
   mHH[16]=MMZ*mHH[5];
   mHH[17]=mHH[3]*MMZ;
   mHH[18]=mHH[12]*mHH[18];
   mHH[15]=2*mHH[18] + 1./2.*mHH[17] + mHH[16] + mHH[9] + mHH[15];
   mHH[15]=mHH[4]*mHH[15];
   mHH[16]=mHH[3]*MMH;
   mHH[13]=mHH[13] + 1./4.*mHH[16];
   mHH[13]=mHH[12]*mHH[13];
   mHH[13]=3*mHH[15] + 1./2.*mHH[13] - mHH[5] - 1./2.*mHH[3];
   mHH[13]=mHH[11]*mHH[13];
   mHH[15]= - mHH[4]*MMZ*mHH[5];

      return mHH[13] + mHH[14] + 3*mHH[15];
}
