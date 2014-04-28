#include <HH.hpp>
std::complex<long double> HH::m11(size_t nG)
{     
      
      
    std::complex<long double> mHH[16];

    mHH[1]=pow(CW,-1);
    mHH[2]=pow(MMH,-1);
    mHH[3]=pow(MMZ,-1);
    mHH[4]=pow(SW,-1);
    mHH[5]=Tsil::B(MMt,MMt,MMH,mu2);
    mHH[6]=Tsil::A(MMt,mu2);
    mHH[7]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mHH[8]=Tsil::Aeps(MMt,mu2);
    mHH[9]=prottttt0->M(0);
    mHH[10]=prottttt0->Vxzuv(0);
    mHH[11]=prottttt0->Suxv(0);
   mHH[12]=2*mHH[10];
   mHH[13]=mHH[12] + mHH[9];
   mHH[14]=MMt*mHH[2];
   mHH[15]=8*mHH[14];
   mHH[13]=mHH[13]*mHH[15];
   mHH[15]=10 + 7*mHH[5];
   mHH[15]=mHH[15]*mHH[5];
   mHH[15]=mHH[15] + 9 - 4*mHH[7];
   mHH[15]=mHH[15]*mHH[2];
   mHH[12]=mHH[12] + 3*mHH[9];
   mHH[12]=mHH[15] - mHH[13] + 2*mHH[12];
   mHH[13]=2*MMt;
   mHH[12]=mHH[12]*mHH[13];
   mHH[13]=mHH[9]*MMH;
   mHH[15]=mHH[2]*mHH[11];
   mHH[13]=mHH[13] - mHH[15];
   mHH[15]=4 + 3*mHH[5];
   mHH[15]=mHH[15]*mHH[5];
   mHH[16]=mHH[8]*mHH[2];
   mHH[12]=mHH[12] - mHH[15] - 12 - 12*mHH[16] - 2*mHH[13];
   mHH[13]=pow(mHH[1],2);
   mHH[15]=pow(mHH[4],2);
   mHH[13]=mHH[13] + mHH[15];
   mHH[13]=mHH[13]*mHH[3];
   mHH[12]= - mHH[12]*MMt*mHH[13];
   mHH[14]= - 1 + 6*mHH[14];
   mHH[14]= - mHH[14]*mHH[5];
   mHH[15]= - mHH[6]*mHH[2];
   mHH[14]=4*mHH[15] + mHH[14];
   mHH[13]=mHH[6]*mHH[13]*mHH[14];
   mHH[12]=6*mHH[13] + mHH[12];

      return 2*mHH[12];
}
