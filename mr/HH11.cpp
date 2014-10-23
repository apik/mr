#include <HH.hpp>
std::complex<long double> HH::m11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHH[18], mHHret;

    mHH[1]=double(nH);
    mHH[2]=pow(CW,-1);
    mHH[3]=pow(MMH,-1);
    mHH[4]=pow(MMZ,-1);
    mHH[5]=pow(SW,-1);
    mHH[6]=Tsil::B(MMt,MMt,MMH,mu2);
    mHH[7]=Tsil::A(MMt,mu2);
    mHH[8]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mHH[9]=Tsil::Aeps(MMt,mu2);
    mHH[10]=prottttt0->M(0);
    mHH[11]=prottttt0->Vxzuv(0);
    mHH[12]=prottttt0->Suxv(0);
   mHH[13]=2*mHH[11];
   mHH[14]=mHH[13] + mHH[10];
   mHH[15]=8*MMt;
   mHH[14]=mHH[14]*mHH[15];
   mHH[15]=10 + 7*mHH[6];
   mHH[15]=mHH[15]*mHH[6];
   mHH[14]= - mHH[15] + mHH[14] - 9 + 4*mHH[8];
   mHH[14]=mHH[14]*MMt;
   mHH[15]=mHH[7]*mHH[6];
   mHH[14]=mHH[14] - mHH[12] - 18*mHH[15];
   mHH[14]=mHH[14]*MMt;
   mHH[16]=pow(mHH[7],2);
   mHH[14]=mHH[14] - 12*mHH[16];
   mHH[14]=mHH[14]*mHH[4];
   mHH[16]=mHH[4]*MMt;
   mHH[17]=mHH[16]*mHH[9];
   mHH[14]=mHH[14] + 6*mHH[17];
   mHH[17]=2*mHH[3];
   mHH[14]=mHH[14]*mHH[17];
   mHH[13]=mHH[13] + 3*mHH[10];
   mHH[17]=4*MMt;
   mHH[13]=mHH[13]*mHH[17];
   mHH[17]=4 + 3*mHH[6];
   mHH[17]=mHH[17]*mHH[6];
   mHH[13]= - mHH[13] + mHH[17] + 12;
   mHH[13]=mHH[13]*MMt;
   mHH[13]=mHH[13] + 6*mHH[15];
   mHH[13]=mHH[13]*mHH[4];
   mHH[15]=mHH[16]*MMH*mHH[10];
   mHH[13]=mHH[14] + mHH[13] + 2*mHH[15];
   mHH[14]=pow(mHH[5],2);
   mHH[15]=pow(mHH[2],2);
   mHH[14]=mHH[14] + mHH[15];

      mHHret = 2*mHH[14]*mHH[13]*mHH[1];
      return mHHret;
}
