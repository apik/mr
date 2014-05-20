#include <WW.hpp>
std::complex<long double> WW::m11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mWW[23];

    mWW[1]=pow(CW,-1);
    mWW[2]=pow(MMH,-1);
    mWW[3]=pow(MMZ,-1);
    mWW[4]=pow(SW,-1);
    mWW[5]=double(nL + nH);
    mWW[6]=Tsil::B(0,0,MMW,mu2);
    mWW[7]=prot00000->M(0);
    mWW[8]=Tsil::B(0,MMt,MMW,mu2);
    mWW[9]=Tsil::A(MMt,mu2);
    mWW[10]=pow(MMt,-1);
    mWW[11]=Tsil::Aeps(MMt,mu2);
    mWW[12]=prot00tt0->M(0);
    mWW[13]=prot00tt0->Tuxv(0);
   mWW[14]=mWW[12] - mWW[7];
   mWW[15]=2./3.*MMZ;
   mWW[14]=mWW[14]*mWW[15];
   mWW[15]=4./3.*mWW[13];
   mWW[16]=mWW[10]*mWW[9];
   mWW[17]=1 - 2*mWW[16];
   mWW[18]= - mWW[12]*MMt;
   mWW[19]=mWW[11]*mWW[10];
   mWW[16]=2./3.*mWW[8] + 1 + 4./3.*mWW[16];
   mWW[16]=mWW[8]*mWW[16];
   mWW[16]=mWW[14] + mWW[15] + mWW[16] + 4./3.*mWW[19] + mWW[18] + 2./3.
   *mWW[17] - mWW[6];
   mWW[17]=pow(mWW[4],2);
   mWW[16]=mWW[16]*mWW[17];
   mWW[14]= - mWW[14] + mWW[16];
   mWW[16]= - 1 + 1./3.*mWW[8];
   mWW[16]=mWW[16]*MMt;
   mWW[16]=mWW[16] + 2./3.*mWW[9];
   mWW[18]=10*mWW[8];
   mWW[16]=mWW[16]*mWW[18];
   mWW[18]=pow(mWW[9],2);
   mWW[19]=pow(MMt,2);
   mWW[20]= - 3*mWW[18] + 4*mWW[19];
   mWW[21]=16*mWW[2];
   mWW[20]=mWW[20]*mWW[21];
   mWW[21]=mWW[13]*MMt;
   mWW[21]=mWW[21] + mWW[11];
   mWW[22]=13./2.*MMt + 4*mWW[9];
   mWW[23]=mWW[18]*mWW[10];
   mWW[16]=mWW[16] - 5./3.*mWW[22] + mWW[20] + 8*mWW[23] + 8./3.*
   mWW[21];
   mWW[20]=pow(mWW[1],2);
   mWW[21]=mWW[20] + mWW[17];
   mWW[16]=mWW[16]*mWW[21];
   mWW[21]=mWW[12]*pow(MMt,3);
   mWW[18]=mWW[21] + mWW[18];
   mWW[21]= - 5*mWW[9] - 4*mWW[11];
   mWW[21]=MMt*mWW[21];
   mWW[18]=2*mWW[18] + mWW[21];
   mWW[15]=mWW[15] + 3*mWW[8];
   mWW[15]=mWW[19]*mWW[15];
   mWW[15]= - mWW[15] + 1./3.*mWW[18];
   mWW[18]=mWW[20] + 1;
   mWW[18]=mWW[20]*mWW[18];
   mWW[18]=mWW[18] + mWW[17];
   mWW[15]=mWW[3]*mWW[15]*mWW[18];
   mWW[15]=2*mWW[15] + mWW[16];
   mWW[15]=mWW[3]*mWW[15];
   mWW[14]=4*mWW[14] + mWW[15];
   mWW[14]=double(nH)*mWW[14];
   mWW[15]=MMZ*mWW[7];
   mWW[15]=8./3.*mWW[15];
   mWW[16]=mWW[15] + 31./3. + 4*mWW[6];
   mWW[16]=mWW[16]*mWW[17];
   mWW[15]= - mWW[15] + mWW[16];
   mWW[15]=mWW[5]*mWW[15];

      return mWW[14] + mWW[15];
}
