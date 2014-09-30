#include <tt.hpp>
std::complex<long double> tt::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myt[21];

    myt[1]=pow(CW,-1);
    myt[2]=pow(MMZ,-1);
    myt[3]=pow(SW,-1);
    myt[4]=Tsil::B(MMH,MMt,MMt,mu2);
    myt[5]=Tsil::B(MMZ,MMt,MMt,mu2);
    myt[6]=pow(MMt,-1);
    myt[7]=Tsil::B(0,MMW,MMt,mu2);
    myt[8]=Tsil::A(MMH,mu2);
    myt[9]=Tsil::A(MMZ,mu2);
    myt[10]=Tsil::A(MMW,mu2);
    myt[11]=Tsil::A(MMt,mu2);
    myt[12]=1/( - MMW + MMH);
   myt[13]=1./4.*myt[7];
   myt[14]= - myt[13] + 3./4.;
   myt[14]=MMt*myt[14];
   myt[15]= - myt[8] + 3*myt[9] + 1./2.*MMH;
   myt[16]= - MMt + 1./4.*MMH;
   myt[16]=myt[16]*myt[4];
   myt[14]=myt[16] - 7./4.*myt[10] + myt[11] - 1./4.*myt[15] + myt[14];
   myt[15]=pow(myt[3],2);
   myt[16]= - myt[9] + myt[10];
   myt[16]=myt[16]*myt[15];
   myt[16]=3./4.*myt[16] - myt[14];
   myt[16]=myt[2]*myt[16];
   myt[17]= - myt[12]*myt[10];
   myt[18]=myt[8]*myt[12];
   myt[17]=myt[18] - 1./2. + myt[17];
   myt[17]=myt[5] + myt[7] + 3*myt[17];
   myt[18]=myt[5]*MMZ;
   myt[19]= - myt[18] + myt[11] - myt[9];
   myt[20]=myt[7]*MMZ;
   myt[20]= - myt[20] - myt[10] + 1./2.*myt[19];
   myt[20]=myt[6]*myt[20];
   myt[16]=1./2.*myt[20] + myt[16] + 1./4.*myt[17];
   myt[15]=myt[16]*myt[15];
   myt[14]= - myt[2]*myt[14];
   myt[14]= - 7./36.*myt[5] + 7./72. + myt[14];
   myt[16]=pow(myt[1],2);
   myt[14]=myt[14]*myt[16];
   myt[14]=myt[14] + myt[15];
   myt[15]=myt[19]*myt[16];
   myt[16]=myt[18] + myt[9] + 2*myt[11];
   myt[13]=MMZ*myt[13];
   myt[13]=17./72.*myt[15] + myt[13] + 4./9.*myt[16];
   myt[13]=myt[6]*myt[13];
   myt[15]= - 1 + myt[5];

      return myt[13] + 1./2.*myt[14] + 8./9.*myt[15];
}
