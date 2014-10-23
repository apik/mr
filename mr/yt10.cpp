#include <tt.hpp>
std::complex<long double> tt::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myt[21], mytret;

    myt[1]=pow(CW,-1);
    myt[2]=pow(MMZ,-1);
    myt[3]=pow(SW,-1);
    myt[4]=double(nH);
    myt[5]=Tsil::A(MMt,mu2);
    myt[6]=Tsil::B(MMH,MMt,MMt,mu2);
    myt[7]=Tsil::B(MMZ,MMt,MMt,mu2);
    myt[8]=pow(MMt,-1);
    myt[9]=Tsil::A(MMH,mu2);
    myt[10]=Tsil::A(MMZ,mu2);
    myt[11]=Tsil::A(MMW,mu2);
    myt[12]=std::real(Tsil::B(0,MMW,MMt,mu2));
    myt[13]=1/( - MMW + MMH);
   myt[14]=MMZ*myt[7];
   myt[14]=myt[14] + myt[10];
   myt[15]=myt[14] - myt[5];
   myt[16]=myt[12]*MMZ;
   myt[16]= - myt[16] - myt[11] - 1./2.*myt[15];
   myt[17]=pow(myt[3],2);
   myt[16]=myt[16]*myt[17];
   myt[18]=pow(myt[1],2);
   myt[15]= - myt[15]*myt[18];
   myt[14]=2*myt[5] + myt[14];
   myt[19]=1./4.*myt[12];
   myt[20]=MMZ*myt[19];
   myt[14]=17./72.*myt[15] + 1./4.*myt[16] + 4./9.*myt[14] + myt[20];
   myt[14]=myt[8]*myt[14];
   myt[15]=myt[6] - 1./2.;
   myt[16]=1./2.*MMH;
   myt[15]=myt[15]*myt[16];
   myt[15]=myt[5] - myt[15] + 3./2.*myt[10] + 7./2.*myt[11];
   myt[16]=myt[5] + 1./2.*MMt;
   myt[20]=3./2.*myt[4];
   myt[16]=myt[16]*myt[20];
   myt[19]=myt[19] + myt[6];
   myt[19]=MMt*myt[19];
   myt[15]= - myt[16] - 1./4.*myt[9] + 1./2.*myt[15] + myt[19];
   myt[16]=myt[11] - myt[10];
   myt[16]=myt[16]*myt[17];
   myt[16]=3./4.*myt[16] + myt[15];
   myt[16]=myt[2]*myt[16];
   myt[19]=myt[9] - myt[11];
   myt[19]=3*myt[19];
   myt[19]=myt[13]*myt[19];
   myt[19]=myt[12] - 3./2. + myt[7] + myt[19];
   myt[16]=1./2.*myt[16] + 1./8.*myt[19];
   myt[16]=myt[17]*myt[16];
   myt[15]=myt[2]*myt[15];
   myt[17]=1./2. - myt[7];
   myt[15]=7./36.*myt[17] + myt[15];
   myt[15]=myt[15]*myt[18];
   myt[17]= - 1 + myt[7];

      mytret = myt[14] + 1./2.*myt[15] + myt[16] + 8./9.*myt[17];
      return mytret;
}
