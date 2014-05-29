#include <tt.hpp>
std::complex<long double> tt::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myt[27];

    myt[1]=pow(CW,-1);
    myt[2]=pow(MMZ,-1);
    myt[3]=pow(SW,-1);
    myt[4]=double(nH);
    myt[5]=Tsil::A(MMt,mu2);
    myt[6]=Tsil::B(MMH,MMt,MMt,mu2);
    myt[7]=Tsil::B(MMZ,MMt,MMt,mu2);
    myt[8]=pow(MMt,-1);
    myt[9]=Tsil::B(MMW,MMb,MMt,mu2);
    myt[10]=Tsil::A(MMH,mu2);
    myt[11]=Tsil::A(MMZ,mu2);
    myt[12]=Tsil::A(MMW,mu2);
    myt[13]=Tsil::A(MMb,mu2);
    myt[14]=1/( - MMb + MMt);
    myt[15]=1/( - MMW + MMH);
   myt[16]=1./2.*MMt;
   myt[17]=myt[16] + myt[5];
   myt[18]= - myt[13] + myt[17] - 1./2.*MMb;
   myt[18]=MMb*myt[18]*myt[14];
   myt[17]=myt[18] + myt[17];
   myt[18]=3*myt[4];
   myt[17]=myt[17]*myt[18];
   myt[18]=myt[12] - myt[13];
   myt[19]=1./2.*myt[8];
   myt[20]=myt[18]*myt[19];
   myt[21]=myt[19]*MMb*myt[9];
   myt[20]= - myt[20] + myt[21] - myt[9];
   myt[20]=myt[20]*MMb;
   myt[22]=myt[13] - myt[10];
   myt[16]=myt[16]*myt[9];
   myt[23]=1./4.*MMH;
   myt[16]= - 3./2.*myt[11] + myt[17] - myt[20] - myt[16] - 7./2.*
   myt[12] - myt[23] - myt[5] - 1./2.*myt[22];
   myt[17]=myt[23] - MMt;
   myt[17]=myt[17]*myt[6];
   myt[16]=myt[17] + 1./2.*myt[16];
   myt[17]=pow(myt[1],2);
   myt[20]= - myt[2]*myt[17]*myt[16];
   myt[22]=pow(myt[3],2);
   myt[23]=myt[12] - myt[11];
   myt[23]=myt[23]*myt[22];
   myt[16]=3./4.*myt[23] - myt[16];
   myt[16]=myt[2]*myt[16];
   myt[23]=MMZ*myt[7];
   myt[24]=myt[23] - myt[5];
   myt[25]=myt[9]*MMZ;
   myt[18]= - myt[25] - 1./2.*myt[24] - myt[18];
   myt[18]=myt[8]*myt[18];
   myt[26]=myt[9] - 3./2. + myt[7];
   myt[19]= - myt[11]*myt[19];
   myt[18]=myt[19] + myt[21] + 1./2.*myt[26] + myt[18];
   myt[19]=myt[10] - myt[12];
   myt[19]=myt[15]*myt[19];
   myt[16]=3./4.*myt[19] + 1./2.*myt[18] + myt[16];
   myt[16]=myt[16]*myt[22];
   myt[16]=myt[20] + myt[16];
   myt[18]=2*myt[5] + myt[23];
   myt[19]=4 - 17./8.*myt[17];
   myt[19]=myt[11]*myt[19];
   myt[18]=1./9.*myt[19] + 4./9.*myt[18] + 1./4.*myt[25];
   myt[18]=myt[8]*myt[18];
   myt[19]=1./2. - myt[7];
   myt[20]=myt[8]*myt[24];
   myt[19]=7*myt[19] - 17*myt[20];
   myt[17]=myt[19]*myt[17];
   myt[19]= - 1 + myt[7];

      return 1./2.*myt[16] + 1./72.*myt[17] + myt[18] + 8./9.*myt[19];
}
