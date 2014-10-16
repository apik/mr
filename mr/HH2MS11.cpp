#include <HH.hpp>
std::complex<long double> hh::m11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHH2MS[18], mHH2MSret;

    mHH2MS[1]=double(nH);
    mHH2MS[2]=pow(c,-1);
    mHH2MS[3]=pow(mmH,-1);
    mHH2MS[4]=pow(mmZ,-1);
    mHH2MS[5]=pow(s,-1);
    mHH2MS[6]=Tsil::B(mmt,mmt,mmH,mu2);
    mHH2MS[7]=Tsil::A(mmt,mu2);
    mHH2MS[8]=Tsil::Beps(mmt,mmt,mmH,mu2);
    mHH2MS[9]=Tsil::Aeps(mmt,mu2);
    mHH2MS[10]=prottttt0->M(0);
    mHH2MS[11]=prottttt0->Vxzuv(0);
    mHH2MS[12]=prottttt0->Suxv(0);
   mHH2MS[13]=2*mHH2MS[3];
   mHH2MS[14]=mHH2MS[13]*mHH2MS[11];
   mHH2MS[15]=mHH2MS[10]*mHH2MS[3];
   mHH2MS[14]=mHH2MS[14] + mHH2MS[15];
   mHH2MS[15]=8*mmt;
   mHH2MS[14]=mHH2MS[14]*mHH2MS[15];
   mHH2MS[15]=mHH2MS[8]*mHH2MS[3];
   mHH2MS[15]=mHH2MS[15] - mHH2MS[11];
   mHH2MS[16]=20 + 7*mHH2MS[6];
   mHH2MS[16]=mHH2MS[16]*mHH2MS[6];
   mHH2MS[16]=mHH2MS[16] + 11;
   mHH2MS[16]=mHH2MS[16]*mHH2MS[3];
   mHH2MS[14]=mHH2MS[16] - mHH2MS[14] + 6*mHH2MS[10] - 4*mHH2MS[15];
   mHH2MS[15]=2*mmt;
   mHH2MS[14]=mHH2MS[14]*mHH2MS[15];
   mHH2MS[15]= - mHH2MS[12] + 6*mHH2MS[9];
   mHH2MS[15]=mHH2MS[15]*mHH2MS[13];
   mHH2MS[16]=mHH2MS[6] + 2;
   mHH2MS[16]=mHH2MS[16]*mHH2MS[6];
   mHH2MS[16]=mHH2MS[16] + 4;
   mHH2MS[17]=mHH2MS[10]*mmH;
   mHH2MS[14]=mHH2MS[14] - 3*mHH2MS[16] - mHH2MS[15] - 2*mHH2MS[17];
   mHH2MS[14]=mHH2MS[14]*mmt;
   mHH2MS[15]= - 1 + 3*mHH2MS[6];
   mHH2MS[13]=mHH2MS[15]*mHH2MS[13]*mmt;
   mHH2MS[15]=mHH2MS[7]*mHH2MS[3];
   mHH2MS[13]=mHH2MS[13] + 9*mHH2MS[15];
   mHH2MS[15]=4*mHH2MS[7];
   mHH2MS[13]=mHH2MS[13]*mHH2MS[15];
   mHH2MS[13]=mHH2MS[14] - mHH2MS[13];
   mHH2MS[14]=pow(mHH2MS[2],2);
   mHH2MS[15]=pow(mHH2MS[5],2);
   mHH2MS[14]=mHH2MS[14] + mHH2MS[15];

      mHH2MSret = 2*mHH2MS[14]*mHH2MS[13]*mHH2MS[4]*mHH2MS[1];
      return mHH2MSret;
}
