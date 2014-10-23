#include <HH.hpp>
std::complex<long double> HH::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHH[22], mHHret;

    mHH[1]=double(nH);
    mHH[2]=Tsil::B(MMt,MMt,MMH,mu2);
    mHH[3]=pow(CW,-1);
    mHH[4]=pow(MMH,-1);
    mHH[5]=pow(MMZ,-1);
    mHH[6]=pow(SW,-1);
    mHH[7]=Tsil::A(MMt,mu2);
    mHH[8]=Tsil::B(MMH,MMH,MMH,mu2);
    mHH[9]=Tsil::B(MMZ,MMZ,MMH,mu2);
    mHH[10]=Tsil::B(MMW,MMW,MMH,mu2);
    mHH[11]=Tsil::A(MMH,mu2);
    mHH[12]=Tsil::A(MMZ,mu2);
    mHH[13]=Tsil::A(MMW,mu2);
   mHH[14]=pow(mHH[6],2);
   mHH[15]=pow(mHH[3],2);
   mHH[16]=mHH[14] + mHH[15];
   mHH[17]=MMt*mHH[2];
   mHH[17]=mHH[17] + mHH[7];
   mHH[18]=mHH[1]*MMt;
   mHH[17]= - mHH[5]*mHH[18]*mHH[17]*mHH[16];
   mHH[19]=1./2.*mHH[9];
   mHH[15]=mHH[19]*mHH[15];
   mHH[20]=mHH[15] - mHH[10];
   mHH[20]=MMZ*mHH[20];
   mHH[19]=mHH[19] + mHH[10];
   mHH[21]=MMZ*mHH[19];
   mHH[21]=mHH[21] + mHH[13];
   mHH[21]=mHH[21]*mHH[14];
   mHH[17]=2*mHH[17] + mHH[21] + mHH[20];
   mHH[20]=3*mHH[4];
   mHH[17]=mHH[17]*mHH[20];
   mHH[21]=mHH[19] + 9./2.*mHH[8];
   mHH[21]=mHH[21]*MMH;
   mHH[21]=mHH[21] + 3*mHH[11];
   mHH[18]=mHH[18]*mHH[2];
   mHH[18]=mHH[13] + 3*mHH[18] + 1./2.*mHH[21];
   mHH[21]=1./2.*mHH[5];
   mHH[18]=mHH[21]*mHH[18];
   mHH[20]=mHH[21] + mHH[20];
   mHH[20]=mHH[12]*mHH[20];
   mHH[18]=1./2.*mHH[20] + mHH[18];
   mHH[16]=mHH[16]*mHH[18];
   mHH[14]= - mHH[19]*mHH[14];

      mHHret = mHH[14] - mHH[15] + mHH[16] + mHH[17];
      return mHHret;
}
