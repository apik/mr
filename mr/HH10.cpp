#include <HH.hpp>
std::complex<long double> HH::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHH[23];

    mHH[1]=double(nH);
    mHH[2]=Tsil::B(MMt,MMt,MMH,mu2);
    mHH[3]=pow(CW,-1);
    mHH[4]=pow(MMH,-1);
    mHH[5]=pow(MMZ,-1);
    mHH[6]=pow(SW,-1);
    mHH[7]=Tsil::B(MMb,MMb,MMH,mu2);
    mHH[8]=Tsil::A(MMt,mu2);
    mHH[9]=Tsil::A(MMb,mu2);
    mHH[10]=Tsil::B(MMH,MMH,MMH,mu2);
    mHH[11]=Tsil::B(MMZ,MMZ,MMH,mu2);
    mHH[12]=Tsil::B(MMW,MMW,MMH,mu2);
    mHH[13]=Tsil::A(MMH,mu2);
    mHH[14]=Tsil::A(MMZ,mu2);
    mHH[15]=Tsil::A(MMW,mu2);
   mHH[16]=MMb*mHH[4];
   mHH[17]=mHH[16]*mHH[9];
   mHH[18]=MMt*mHH[4];
   mHH[19]=mHH[18]*mHH[8];
   mHH[17]=mHH[17] + mHH[19];
   mHH[16]= - 1./2. + 2*mHH[16];
   mHH[16]=mHH[16]*mHH[7]*MMb;
   mHH[18]= - 1./2. + 2*mHH[18];
   mHH[18]=mHH[18]*MMt*mHH[2];
   mHH[16]=mHH[18] + mHH[16] + 2*mHH[17];
   mHH[17]=3*mHH[1];
   mHH[16]=mHH[17]*mHH[16];
   mHH[17]=mHH[14] + 3*mHH[13];
   mHH[16]= - 1./4.*mHH[17] + mHH[16];
   mHH[16]=mHH[16]*mHH[5];
   mHH[17]=1./2.*mHH[11];
   mHH[18]=mHH[17] + mHH[12];
   mHH[19]=9./2.*mHH[10] + mHH[18];
   mHH[20]=1./4.*MMH;
   mHH[19]=mHH[20]*mHH[19]*mHH[5];
   mHH[20]=mHH[14]*mHH[4];
   mHH[16]= - mHH[19] + mHH[16] - 3./2.*mHH[20];
   mHH[19]=3*mHH[4];
   mHH[20]=mHH[19]*MMZ;
   mHH[21]=mHH[20] - 1;
   mHH[18]=mHH[21]*mHH[18];
   mHH[22]=1./2.*mHH[5];
   mHH[19]=mHH[19] + mHH[22];
   mHH[19]=mHH[15]*mHH[19];
   mHH[18]=mHH[19] + mHH[18] - mHH[16];
   mHH[18]=mHH[18]*pow(mHH[6],2);
   mHH[17]=mHH[21]*mHH[17];
   mHH[19]=mHH[15]*mHH[22];
   mHH[16]=mHH[17] + mHH[19] - mHH[16];
   mHH[16]=mHH[16]*pow(mHH[3],2);
   mHH[17]= - mHH[12]*mHH[20];

      return mHH[16] + mHH[17] + mHH[18];
}
