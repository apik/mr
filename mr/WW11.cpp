#include <WW.hpp>
std::complex<long double> WW::m11(size_t nG = 3)
{     
      
      
    std::complex<long double> mWW[24];

    mWW[1]=pow(CW,-1);
    mWW[2]=pow(MMH,-1);
    mWW[3]=pow(MMZ,-1);
    mWW[4]=pow(SW,-1);
    mWW[5]=double(nG);
    mWW[6]=Tsil::B(0,0,MMW,mu2);
    mWW[7]=prot00000->M(0);
    mWW[8]=Tsil::B(0,MMt,MMW,mu2);
    mWW[9]=Tsil::A(MMt,mu2);
    mWW[10]=pow(MMt,-1);
    mWW[11]=Tsil::Aeps(MMt,mu2);
    mWW[12]=prot00tt0->M(0);
    mWW[13]=prot00tt0->Tuxv(0);
   mWW[14]=2*mWW[10];
   mWW[15]=mWW[11] - mWW[9];
   mWW[15]=mWW[15]*mWW[14];
   mWW[16]=2*mWW[13];
   mWW[15]=mWW[15] + 1 + mWW[16];
   mWW[17]=MMt*mWW[12];
   mWW[18]=mWW[10]*mWW[9];
   mWW[18]=2./3.*mWW[8] + 1 + 4./3.*mWW[18];
   mWW[18]=mWW[8]*mWW[18];
   mWW[19]=mWW[5] - 1;
   mWW[20]=mWW[6]*mWW[19];
   mWW[15]=mWW[20] + mWW[18] + 2./3.*mWW[15] - mWW[17];
   mWW[18]=mWW[19]*mWW[7];
   mWW[18]=mWW[18] + mWW[12];
   mWW[18]=MMZ*mWW[18];
   mWW[18]=8./3.*mWW[18];
   mWW[15]=31./3.*mWW[5] + 4*mWW[15] + mWW[18];
   mWW[19]=pow(mWW[4],2);
   mWW[15]=mWW[19]*mWW[15];
   mWW[20]= - 65./2. + 8*mWW[13];
   mWW[21]=MMt*mWW[2];
   mWW[20]=64*mWW[21] + 1./3.*mWW[20];
   mWW[20]=mWW[20]*MMt;
   mWW[14]=mWW[14] - 12*mWW[2];
   mWW[21]=pow(mWW[9],2);
   mWW[14]=mWW[21]*mWW[14];
   mWW[22]=5*mWW[9];
   mWW[23]= - mWW[22] + 2*mWW[11];
   mWW[14]=1./3.*mWW[23] + mWW[14];
   mWW[23]= - 1 + 1./3.*mWW[8];
   mWW[23]=MMt*mWW[23];
   mWW[23]=2./3.*mWW[9] + mWW[23];
   mWW[24]=10*mWW[8];
   mWW[23]=mWW[23]*mWW[24];
   mWW[14]=mWW[23] + mWW[20] + 4*mWW[14];
   mWW[20]=pow(mWW[1],2);
   mWW[23]=mWW[20] + mWW[19];
   mWW[14]=mWW[14]*mWW[23];
   mWW[23]= - mWW[20] - 1;
   mWW[20]=mWW[20]*mWW[23];
   mWW[19]=mWW[20] - mWW[19];
   mWW[16]=mWW[16] - mWW[17];
   mWW[17]=2*MMt;
   mWW[16]=mWW[16]*mWW[17];
   mWW[16]=mWW[16] + mWW[22] + 4*mWW[11];
   mWW[16]=mWW[16]*MMt;
   mWW[16]=mWW[16] - 2*mWW[21];
   mWW[17]=3*mWW[8];
   mWW[17]=mWW[17]*pow(MMt,2);
   mWW[16]=mWW[17] + 1./3.*mWW[16];
   mWW[16]=mWW[3]*mWW[16]*mWW[19];
   mWW[14]=2*mWW[16] + mWW[14];
   mWW[14]=mWW[3]*mWW[14];

      return mWW[14] + mWW[15] - mWW[18];
}
