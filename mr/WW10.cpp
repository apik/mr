#include <WW.hpp>
std::complex<long double> WW::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mWW[20], mWWret;

    mWW[1]=pow(CW,-1);
    mWW[2]=pow(MMH,-1);
    mWW[3]=pow(MMZ,-1);
    mWW[4]=pow(SW,-1);
    mWW[5]=double(nL + nH);
    mWW[6]=std::real(Tsil::B(0,0,MMW,mu2));
    mWW[7]=Tsil::B(MMW,MMH,MMW,mu2);
    mWW[8]=Tsil::B(MMW,MMZ,MMW,mu2);
    mWW[9]=Tsil::B(0,MMt,MMW,mu2);
    mWW[10]=Tsil::A(MMH,mu2);
    mWW[11]=Tsil::A(MMZ,mu2);
    mWW[12]=Tsil::A(MMW,mu2);
    mWW[13]=Tsil::A(MMt,mu2);
   mWW[14]=MMH*mWW[7];
   mWW[14]=mWW[14] + mWW[10] - mWW[12];
   mWW[15]=1./6.*MMH;
   mWW[14]=mWW[14]*mWW[15];
   mWW[15]=MMt*mWW[9];
   mWW[15]=mWW[15] + mWW[13];
   mWW[15]=mWW[15]*MMt;
   mWW[14]=mWW[14] - mWW[15];
   mWW[15]=1./2.*mWW[3];
   mWW[15]=mWW[14]*mWW[15];
   mWW[16]=mWW[7] + 1./2.;
   mWW[17]=1./3.*MMH;
   mWW[16]=mWW[16]*mWW[17];
   mWW[17]=mWW[2]*mWW[13];
   mWW[17]=6*mWW[17] - 1 + 1./2.*mWW[9];
   mWW[17]=mWW[17]*MMt;
   mWW[15]= - mWW[16] + mWW[13] + 1./2.*mWW[10] + mWW[15] - mWW[17];
   mWW[16]=3./4.*mWW[11] + 35./12.*mWW[12] + mWW[15];
   mWW[16]=mWW[3]*mWW[16];
   mWW[14]=mWW[3]*mWW[14];
   mWW[17]= - mWW[12] + mWW[11];
   mWW[14]=1./6.*mWW[17] + mWW[14];
   mWW[14]=mWW[3]*mWW[14];
   mWW[14]=1./6.*mWW[8] + mWW[14];
   mWW[17]=pow(mWW[1],2);
   mWW[14]=mWW[14]*mWW[17];
   mWW[18]=mWW[11]*mWW[2];
   mWW[18]=3./2.*mWW[18];
   mWW[19]=mWW[2]*MMZ;
   mWW[14]=1./2.*mWW[14] + mWW[16] + 17./12.*mWW[8] + mWW[18] - 1./6.
    + mWW[19];
   mWW[14]=mWW[14]*mWW[17];
   mWW[15]= - 5./4.*mWW[11] - 13./12.*mWW[12] + mWW[15];
   mWW[15]=mWW[3]*mWW[15];
   mWW[16]=mWW[12]*mWW[2];
   mWW[16]=mWW[19] + mWW[16];
   mWW[17]= - 1 + 4./3.*mWW[5];
   mWW[17]=mWW[6]*mWW[17];
   mWW[15]=mWW[15] - 33./4.*mWW[8] + mWW[18] + mWW[7] + mWW[17] - 4./9.
   *mWW[5] + mWW[9] - 59./18. + 3*mWW[16];
   mWW[15]=mWW[15]*pow(mWW[4],2);
   mWW[16]=2*mWW[8] - 2 - mWW[19];

      mWWret = mWW[14] + mWW[15] + 2*mWW[16];
      return mWWret;
}
