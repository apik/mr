#include <HH.hpp>
std::complex<long double> HH::m11(size_t nG)
{     
      
      
    std::complex<long double> mHH[16];

    mHH[1]=pow(MMH,-1);
    mHH[2]=pow(MMW,-1);
    mHH[3]=Tsil::B(MMt,MMt,MMH,mu2);
    mHH[4]=Tsil::A(MMt,mu2);
    mHH[5]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mHH[6]=Tsil::Aeps(MMt,mu2);
    mHH[7]=prottttt0->M(0);
    mHH[8]=prottttt0->Vxzuv(0);
    mHH[9]=prottttt0->Suxv(0);
    mHH[10]=1/(1 - pow(MMZ,-1)*MMW);
    mHH[11]=pow(MMZ,-1);
   mHH[12]= - 2*mHH[8] - 3*mHH[7];
   mHH[13]=mHH[10]*mHH[11]*mHH[12];
   mHH[12]=mHH[2]*mHH[12];
   mHH[12]=mHH[13] + mHH[12];
   mHH[13]= - 10 - 7*mHH[3];
   mHH[13]=mHH[3]*mHH[13];
   mHH[13]=mHH[13] - 9 + 4*mHH[5];
   mHH[14]=mHH[10]*mHH[11]*mHH[13];
   mHH[13]=mHH[2]*mHH[13];
   mHH[13]=mHH[14] + mHH[13];
   mHH[13]=mHH[1]*mHH[13];
   mHH[14]=2*mHH[8] + mHH[7];
   mHH[15]=mHH[10]*mHH[11]*mHH[14];
   mHH[14]=mHH[2]*mHH[14];
   mHH[14]=mHH[15] + mHH[14];
   mHH[14]=MMt*mHH[1]*mHH[14];
   mHH[12]=8*mHH[14] + 2*mHH[12] + mHH[13];
   mHH[12]=MMt*mHH[12];
   mHH[13]=mHH[7]*MMH;
   mHH[13]=6 + mHH[13];
   mHH[14]=4 + 3*mHH[3];
   mHH[14]=mHH[3]*mHH[14];
   mHH[13]=2*mHH[13] + mHH[14];
   mHH[14]=mHH[10]*mHH[11]*mHH[13];
   mHH[13]=mHH[2]*mHH[13];
   mHH[15]= - mHH[3]*mHH[4];
   mHH[15]=18*mHH[15] - mHH[9] + 6*mHH[6];
   mHH[16]=mHH[10]*mHH[11]*mHH[15];
   mHH[15]=mHH[2]*mHH[15];
   mHH[15]=mHH[16] + mHH[15];
   mHH[15]=mHH[1]*mHH[15];
   mHH[12]=2*mHH[12] + 2*mHH[15] + mHH[14] + mHH[13];
   mHH[12]=MMt*mHH[12];
   mHH[13]=pow(mHH[4],2);
   mHH[14]= - mHH[10]*mHH[11]*mHH[13];
   mHH[13]= - mHH[2]*mHH[13];
   mHH[13]=mHH[14] + mHH[13];
   mHH[13]=mHH[1]*mHH[13];
   mHH[14]=mHH[3]*mHH[4];
   mHH[15]=mHH[10]*mHH[11]*mHH[14];
   mHH[14]=mHH[2]*mHH[14];
   mHH[13]=4*mHH[13] + mHH[15] + mHH[14];
   mHH[12]=6*mHH[13] + mHH[12];

      return 2*mHH[12];
}
