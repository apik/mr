#include <ZZ.hpp>
std::complex<long double> ZZ::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mZZ[24];

    mZZ[1]=pow(CW,-1);
    mZZ[2]=pow(MMH,-1);
    mZZ[3]=pow(MMZ,-1);
    mZZ[4]=pow(SW,-1);
    mZZ[5]=double(nL + nH);
    mZZ[6]=Tsil::B(0,0,MMZ,mu2);
    mZZ[7]=double(nH);
    mZZ[8]=Tsil::B(MMt,MMt,MMZ,mu2);
    mZZ[9]=Tsil::A(MMt,mu2);
    mZZ[10]=Tsil::B(MMZ,MMH,MMZ,mu2);
    mZZ[11]=Tsil::B(MMW,MMW,MMZ,mu2);
    mZZ[12]=Tsil::A(MMH,mu2);
    mZZ[13]=Tsil::A(MMZ,mu2);
    mZZ[14]=Tsil::A(MMW,mu2);
   mZZ[15]=pow(mZZ[1],2);
   mZZ[16]=pow(mZZ[4],2);
   mZZ[17]=17./9.*mZZ[15] + mZZ[16] - 32./9.;
   mZZ[18]=1./2.*mZZ[16];
   mZZ[19]=7./18.*mZZ[15] - 32./9. - mZZ[18];
   mZZ[19]=mZZ[8]*mZZ[19];
   mZZ[20]=mZZ[15] + mZZ[16];
   mZZ[21]=mZZ[2]*mZZ[9]*mZZ[20];
   mZZ[19]= - 6*mZZ[21] + mZZ[19] + mZZ[17];
   mZZ[19]=MMt*mZZ[19];
   mZZ[17]=mZZ[9]*mZZ[17];
   mZZ[17]=mZZ[19] + mZZ[17];
   mZZ[17]=mZZ[7]*mZZ[17];
   mZZ[19]=1./3.*MMH;
   mZZ[21]= - mZZ[19] + mZZ[12];
   mZZ[19]= - mZZ[10]*mZZ[19];
   mZZ[19]=1./2.*mZZ[21] + 1./6.*mZZ[13] + mZZ[19];
   mZZ[19]=mZZ[20]*mZZ[19];
   mZZ[21]= - mZZ[13] + mZZ[12];
   mZZ[21]=MMH*mZZ[20]*mZZ[21];
   mZZ[22]=mZZ[20]*mZZ[10];
   mZZ[23]=pow(MMH,2)*mZZ[22];
   mZZ[21]=mZZ[23] + mZZ[21];
   mZZ[21]=mZZ[3]*mZZ[21];
   mZZ[17]=1./12.*mZZ[21] + mZZ[19] + mZZ[17];
   mZZ[17]=mZZ[3]*mZZ[17];
   mZZ[19]=3*mZZ[16];
   mZZ[21]=mZZ[15] - 2 + mZZ[19];
   mZZ[21]=MMZ*mZZ[21];
   mZZ[20]=mZZ[13]*mZZ[20];
   mZZ[20]=3./2.*mZZ[20] + mZZ[21];
   mZZ[20]=mZZ[2]*mZZ[20];
   mZZ[21]=5./3.*mZZ[15] + mZZ[16] - 8./3.;
   mZZ[23]=4./3.*mZZ[5];
   mZZ[21]=mZZ[21]*mZZ[23];
   mZZ[23]=1./6.*mZZ[15];
   mZZ[24]= - mZZ[21] - mZZ[23] + 8 - 59./6.*mZZ[16];
   mZZ[18]=mZZ[18] - 16./9. + 17./18.*mZZ[15];
   mZZ[18]=mZZ[18]*mZZ[7];
   mZZ[21]=mZZ[21] - mZZ[18];
   mZZ[21]=mZZ[6]*mZZ[21];
   mZZ[15]=1./12.*mZZ[15] + 29./3. - 33./4.*mZZ[16];
   mZZ[15]=mZZ[11]*mZZ[15];
   mZZ[18]=mZZ[8]*mZZ[18];
   mZZ[16]=mZZ[23] + 4 - 5./2.*mZZ[16];
   mZZ[16]=mZZ[3]*mZZ[16];
   mZZ[19]=mZZ[2]*mZZ[19];
   mZZ[16]=mZZ[19] + mZZ[16];
   mZZ[16]=mZZ[14]*mZZ[16];
   mZZ[19]=1 + mZZ[11];
   mZZ[19]=mZZ[19]*pow(CW,2);

      return mZZ[15] + mZZ[16] + mZZ[17] + mZZ[18] + 4*mZZ[19] + 
      mZZ[20] + mZZ[21] + mZZ[22] + 1./3.*mZZ[24];
}
