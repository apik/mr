#include <WW.hpp>
std::complex<long double> WW::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mWW[22];

    mWW[1]=pow(CW,-1);
    mWW[2]=pow(MMH,-1);
    mWW[3]=pow(MMZ,-1);
    mWW[4]=pow(SW,-1);
    mWW[5]=double(nL + nH);
    mWW[6]=Tsil::B(0,0,MMW,mu2);
    mWW[7]=Tsil::B(MMW,MMH,MMW,mu2);
    mWW[8]=Tsil::B(MMW,MMZ,MMW,mu2);
    mWW[9]=Tsil::B(0,MMt,MMW,mu2);
    mWW[10]=Tsil::A(MMH,mu2);
    mWW[11]=Tsil::A(MMZ,mu2);
    mWW[12]=Tsil::A(MMW,mu2);
    mWW[13]=Tsil::A(MMt,mu2);
   mWW[14]= - 1./2.*mWW[9] + 1;
   mWW[14]=MMt*mWW[14];
   mWW[15]=mWW[13]*MMt;
   mWW[16]=mWW[15]*mWW[2];
   mWW[14]= - 6*mWW[16] + mWW[13] + mWW[14];
   mWW[14]=mWW[14]*double(nH);
   mWW[16]=mWW[7] + 1./2.;
   mWW[17]=1./3.*MMH;
   mWW[16]=mWW[16]*mWW[17];
   mWW[14]=mWW[14] - mWW[16];
   mWW[16]=pow(mWW[1],2);
   mWW[17]=1./12.*mWW[16];
   mWW[18]= - mWW[12] + mWW[11];
   mWW[18]=mWW[18]*mWW[17];
   mWW[19]=3./2.*mWW[11];
   mWW[20]=mWW[19] + mWW[10] + 35./6.*mWW[12];
   mWW[18]=mWW[18] + 1./2.*mWW[20] + mWW[14];
   mWW[18]=mWW[16]*mWW[18];
   mWW[20]= - 5./2.*mWW[11] + mWW[10] - 13./6.*mWW[12];
   mWW[14]=1./2.*mWW[20] + mWW[14];
   mWW[20]=pow(mWW[4],2);
   mWW[14]=mWW[14]*mWW[20];
   mWW[21]=mWW[9]*pow(MMt,2);
   mWW[15]=mWW[21] + mWW[15];
   mWW[15]=mWW[15]*double(nH);
   mWW[21]=MMH*mWW[7];
   mWW[21]= - mWW[21] + mWW[12] - mWW[10];
   mWW[21]=mWW[21]*MMH;
   mWW[15]=mWW[15] + 1./6.*mWW[21];
   mWW[21]=mWW[16] + 1;
   mWW[21]= - mWW[16]*mWW[21]*mWW[15];
   mWW[15]= - mWW[15]*mWW[20];
   mWW[15]=mWW[15] + mWW[21];
   mWW[15]=mWW[3]*mWW[15];
   mWW[14]=1./2.*mWW[15] + mWW[14] + mWW[18];
   mWW[14]=mWW[3]*mWW[14];
   mWW[15]=mWW[12] + MMZ;
   mWW[15]=mWW[19] + 3*mWW[15];
   mWW[15]=mWW[2]*mWW[15];
   mWW[18]= - 59./2. - 4*mWW[5];
   mWW[21]=mWW[6]*mWW[5];
   mWW[22]= - mWW[6] + mWW[9];
   mWW[22]=double(nH)*mWW[22];
   mWW[15]=mWW[22] + mWW[7] + 4./3.*mWW[21] - 33./4.*mWW[8] + 1./9.*
   mWW[18] + mWW[15];
   mWW[15]=mWW[15]*mWW[20];
   mWW[17]=mWW[17] + 17./12.;
   mWW[17]=mWW[8]*mWW[17];
   mWW[18]=mWW[2]*MMZ;
   mWW[19]=mWW[2]*mWW[19];
   mWW[17]=mWW[19] - 1./6. + mWW[18] + mWW[17];
   mWW[16]=mWW[17]*mWW[16];
   mWW[17]=2*mWW[8] - 2 - mWW[18];

      return mWW[14] + mWW[15] + mWW[16] + 2*mWW[17];
}
