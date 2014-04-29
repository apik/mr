#include <HH.hpp>
std::complex<long double> HH::m11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHH[17];

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
   mHH[12]=10 + 7*mHH[5];
   mHH[12]=mHH[12]*mHH[5];
   mHH[12]=mHH[12] + 9 - 4*mHH[7];
   mHH[12]=mHH[12]*mHH[2];
   mHH[13]=MMt*mHH[2];
   mHH[14]= - 4 + 16*mHH[13];
   mHH[14]=mHH[14]*mHH[10];
   mHH[12]=mHH[12] - mHH[14];
   mHH[14]=2*MMt;
   mHH[12]=mHH[12]*mHH[14];
   mHH[15]=mHH[6]*mHH[5];
   mHH[16]= - mHH[11] - 18*mHH[15] + 6*mHH[8];
   mHH[17]=2*mHH[2];
   mHH[16]=mHH[16]*mHH[17];
   mHH[17]=4 + 3*mHH[5];
   mHH[17]=mHH[17]*mHH[5];
   mHH[12]=mHH[12] - mHH[16] - mHH[17] - 12;
   mHH[12]=mHH[12]*MMt;
   mHH[16]=4*mHH[2];
   mHH[16]=mHH[16]*pow(mHH[6],2);
   mHH[15]=mHH[16] - mHH[15];
   mHH[13]= - 3 + 4*mHH[13];
   mHH[13]=mHH[13]*mHH[14];
   mHH[13]=mHH[13] + MMH;
   mHH[13]=mHH[9]*MMt*mHH[13];
   mHH[12]=2*mHH[13] - mHH[12] - 6*mHH[15];
   mHH[13]=pow(mHH[1],2);
   mHH[14]=pow(mHH[4],2);
   mHH[13]=mHH[13] + mHH[14];

      return 2*nH*mHH[13]*mHH[12]*mHH[3];
}
