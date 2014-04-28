#include <WW.hpp>
std::complex<long double> WW::m10(size_t nG)
{     
      
      
    std::complex<long double> mWW[20];

    mWW[1]=pow(MMH,-1);
    mWW[2]=pow(MMW,-1);
    mWW[3]=Tsil::B(MMW,MMH,MMW,mu2);
    mWW[4]=pow(MMZ,-1);
    mWW[5]=Tsil::B(MMW,MMZ,MMW,mu2);
    mWW[6]=Tsil::B(0,MMt,MMW,mu2);
    mWW[7]=Tsil::A(MMH,mu2);
    mWW[8]=Tsil::A(MMZ,mu2);
    mWW[9]=Tsil::A(MMW,mu2);
    mWW[10]=Tsil::A(MMt,mu2);
    mWW[11]=1/(1 - pow(MMZ,-1)*MMW);
    mWW[12]=double(nG);
    mWW[13]=Tsil::B(0,0,MMW,mu2);
   mWW[14]=MMZ*mWW[5];
   mWW[14]=mWW[14] + mWW[8] - mWW[9];
   mWW[14]=MMZ*mWW[14];
   mWW[15]= - MMt*mWW[6];
   mWW[15]= - mWW[10] + mWW[15];
   mWW[15]=MMt*mWW[15];
   mWW[16]=MMH*mWW[3];
   mWW[16]=mWW[16] + mWW[7] - mWW[9];
   mWW[16]=1./6.*MMH*mWW[16];
   mWW[14]=mWW[16] + 1./6.*mWW[14] + mWW[15];
   mWW[14]=mWW[2]*mWW[14];
   mWW[17]= - 1 + 17./2.*mWW[5];
   mWW[18]=mWW[1]*mWW[8];
   mWW[17]=1./3.*mWW[17] + 3*mWW[18];
   mWW[18]=MMZ*mWW[1];
   mWW[17]=1./2.*mWW[17] + mWW[18];
   mWW[17]=MMZ*mWW[17];
   mWW[15]=mWW[15] + mWW[16];
   mWW[15]=1./2.*mWW[4]*mWW[15];
   mWW[16]=3./2.*mWW[8] + mWW[7];
   mWW[19]= - mWW[1]*mWW[10];
   mWW[19]=6*mWW[19] + 1 - 1./2.*mWW[6];
   mWW[19]=MMt*mWW[19];
   mWW[20]= - 1./2. - mWW[3];
   mWW[20]=1./3.*MMH*mWW[20];
   mWW[14]=1./2.*mWW[14] + mWW[15] + mWW[20] + mWW[19] + mWW[17] + 35./
   12.*mWW[9] + 1./2.*mWW[16] + mWW[10];
   mWW[14]=mWW[2]*mWW[14];
   mWW[16]= - 5./2.*mWW[8] + mWW[7];
   mWW[15]=mWW[15] + mWW[20] + mWW[19] - 13./12.*mWW[9] + 1./2.*mWW[16]
    + mWW[10];
   mWW[15]=mWW[4]*mWW[15];
   mWW[16]= - 1./3. + mWW[13];
   mWW[16]=mWW[12]*mWW[16];
   mWW[17]=1./2.*mWW[8] + mWW[9];
   mWW[17]=mWW[1]*mWW[17];
   mWW[15]=mWW[15] + 3*mWW[18] + 3*mWW[17] + mWW[3] + mWW[6] - 33./4.*
   mWW[5] + 4./3.*mWW[16] - 59./18. - mWW[13];
   mWW[15]=mWW[11]*mWW[15];
   mWW[16]= - 1 + mWW[5];
   mWW[17]= - MMZ*mWW[1];
   mWW[16]=2*mWW[16] + mWW[17];

      return mWW[14] + mWW[15] + 2*mWW[16];
}
