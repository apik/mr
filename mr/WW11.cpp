#include <WW.hpp>
std::complex<long double> WW::m11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mWW[12], mWWret;

    mWW[1]=pow(CW,-1);
    mWW[2]=pow(MMH,-1);
    mWW[3]=pow(MMZ,-1);
    mWW[4]=pow(SW,-1);
    mWW[5]=Tsil::B(0,MMt,MMW,mu2);
    mWW[6]=Tsil::A(MMt,mu2);
   mWW[7]=3*mWW[6];
   mWW[8]=mWW[7]*MMt;
   mWW[9]=pow(MMt,2);
   mWW[8]=mWW[8] - mWW[9];
   mWW[8]=mWW[8]*mWW[5];
   mWW[7]=mWW[7] - MMt;
   mWW[10]=mWW[7]*mWW[6];
   mWW[8]=mWW[8] + mWW[10];
   mWW[10]=mWW[8]*mWW[3];
   mWW[11]=MMt + 6*mWW[6];
   mWW[11]=mWW[11]*mWW[6];
   mWW[9]=mWW[11] - mWW[9];
   mWW[11]=4*mWW[2];
   mWW[9]=mWW[9]*mWW[11];
   mWW[11]=mWW[5] - 2;
   mWW[7]=mWW[11]*mWW[7];
   mWW[7]=mWW[10] + mWW[9] + mWW[7];
   mWW[7]=mWW[7]*mWW[3];
   mWW[9]= - pow(mWW[4],2)*mWW[7];
   mWW[10]=pow(mWW[1],2);
   mWW[8]= - mWW[8]*pow(mWW[3],2)*mWW[10];
   mWW[7]= - mWW[7] + mWW[8];
   mWW[7]=mWW[7]*mWW[10];
   mWW[7]=mWW[9] + mWW[7];

      mWWret = 4*mWW[7];
      return mWWret;
}
