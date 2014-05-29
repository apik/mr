#include <WW.hpp>
std::complex<long double> WW::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mWW[24];

    mWW[1]=pow(CW,-1);
    mWW[2]=pow(MMH,-1);
    mWW[3]=pow(MMZ,-1);
    mWW[4]=pow(SW,-1);
    mWW[5]=double(nL + nH);
    mWW[6]=Tsil::B(0,0,MMW,mu2);
    mWW[7]=double(nH);
    mWW[8]=Tsil::B(MMt,MMb,MMW,mu2);
    mWW[9]=Tsil::A(MMt,mu2);
    mWW[10]=Tsil::A(MMb,mu2);
    mWW[11]=Tsil::B(MMW,MMH,MMW,mu2);
    mWW[12]=Tsil::B(MMW,MMZ,MMW,mu2);
    mWW[13]=Tsil::A(MMH,mu2);
    mWW[14]=Tsil::A(MMZ,mu2);
    mWW[15]=Tsil::A(MMW,mu2);
   mWW[16]=mWW[11]*pow(MMH,2);
   mWW[17]=mWW[15]*MMH;
   mWW[18]=mWW[13]*MMH;
   mWW[16]=mWW[18] + mWW[16] - mWW[17];
   mWW[17]=1./2.*mWW[9] - 1./2.*mWW[10];
   mWW[18]=MMt - MMb;
   mWW[17]=mWW[18]*mWW[17];
   mWW[18]= - MMb + 1./2.*MMt;
   mWW[18]=mWW[18]*MMt;
   mWW[19]=pow(MMb,2);
   mWW[18]=mWW[18] + 1./2.*mWW[19];
   mWW[18]=mWW[18]*mWW[8];
   mWW[17]=mWW[18] + mWW[17];
   mWW[17]=mWW[17]*mWW[7];
   mWW[16]=mWW[17] - 1./12.*mWW[16];
   mWW[16]=mWW[16]*mWW[3];
   mWW[17]=MMt + MMb;
   mWW[18]=1./2.*mWW[8];
   mWW[17]=mWW[17]*mWW[18];
   mWW[18]=6*mWW[2];
   mWW[19]=mWW[18]*MMb;
   mWW[19]=mWW[19] - 1;
   mWW[19]=mWW[19]*mWW[10];
   mWW[18]=mWW[18]*mWW[9];
   mWW[18]=mWW[18] - 1;
   mWW[18]=mWW[18]*MMt;
   mWW[17]= - mWW[17] - mWW[19] - mWW[18] + mWW[9] + MMb;
   mWW[17]=mWW[17]*mWW[7];
   mWW[18]=mWW[11] + 1./2.;
   mWW[19]=1./3.*MMH;
   mWW[18]=mWW[18]*mWW[19];
   mWW[17]=mWW[16] - mWW[17] + mWW[18] - 1./2.*mWW[13];
   mWW[18]= - 13./12.*mWW[15] - 5./4.*mWW[14] - mWW[17];
   mWW[18]=mWW[3]*mWW[18];
   mWW[19]=MMZ + mWW[15];
   mWW[20]=3*mWW[2];
   mWW[19]=mWW[20]*mWW[19];
   mWW[20]=mWW[20]*mWW[14];
   mWW[21]= - 59./9. + mWW[20];
   mWW[22]=mWW[7]*mWW[8];
   mWW[23]=4./3.*mWW[5] - mWW[7];
   mWW[23]=mWW[6]*mWW[23];
   mWW[18]=mWW[18] + mWW[23] + mWW[22] - 4./9.*mWW[5] + 1./2.*mWW[21]
    + mWW[11] + mWW[19];
   mWW[19]=pow(mWW[4],2);
   mWW[18]=mWW[18]*mWW[19];
   mWW[17]=35./12.*mWW[15] + 3./4.*mWW[14] - mWW[17];
   mWW[17]=mWW[3]*mWW[17];
   mWW[21]=mWW[14] - mWW[15];
   mWW[16]=1./12.*mWW[21] - mWW[16];
   mWW[16]=mWW[3]*mWW[16];
   mWW[16]=mWW[16] + 1./12.*mWW[12];
   mWW[21]=pow(mWW[1],2);
   mWW[16]=mWW[16]*mWW[21];
   mWW[22]=MMZ*mWW[2];
   mWW[20]= - 1./3. + mWW[20];
   mWW[16]=mWW[16] + 17./12.*mWW[12] + mWW[17] + 1./2.*mWW[20] + 
   mWW[22];
   mWW[16]=mWW[16]*mWW[21];
   mWW[17]= - 2 - mWW[22];
   mWW[19]=4 - 33./4.*mWW[19];
   mWW[19]=mWW[12]*mWW[19];

      return mWW[16] + 2*mWW[17] + mWW[18] + mWW[19];
}
