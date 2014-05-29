#include <tt.hpp>
std::complex<long double> tt::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mtt[25];

    mtt[1]=pow(CW,-1);
    mtt[2]=pow(MMH,-1);
    mtt[3]=pow(SW,-1);
    mtt[4]=double(nH);
    mtt[5]=Tsil::A(MMt,mu2);
    mtt[6]=pow(MMZ,-1);
    mtt[7]=Tsil::A(MMb,mu2);
    mtt[8]=Tsil::B(MMH,MMt,MMt,mu2);
    mtt[9]=Tsil::B(MMZ,MMt,MMt,mu2);
    mtt[10]=pow(MMt,-1);
    mtt[11]=Tsil::B(MMW,MMb,MMt,mu2);
    mtt[12]=Tsil::A(MMH,mu2);
    mtt[13]=Tsil::A(MMZ,mu2);
    mtt[14]=Tsil::A(MMW,mu2);
   mtt[15]=mtt[9]*MMZ;
   mtt[16]=mtt[15] - mtt[5];
   mtt[17]=MMb*mtt[11];
   mtt[18]=mtt[16] - mtt[17];
   mtt[19]=mtt[7] - mtt[14];
   mtt[20]= - mtt[11]*MMZ;
   mtt[18]=mtt[20] + mtt[19] - 1./2.*mtt[18];
   mtt[20]=pow(mtt[3],2);
   mtt[18]=mtt[18]*mtt[20];
   mtt[21]=pow(mtt[1],2);
   mtt[22]=17./72.*mtt[21];
   mtt[16]= - mtt[16]*mtt[22];
   mtt[23]=1./4.*mtt[11];
   mtt[24]=MMZ*mtt[23];
   mtt[15]=mtt[16] + 1./4.*mtt[18] + 8./9.*mtt[5] + mtt[24] + 4./9.*
   mtt[15];
   mtt[15]=mtt[10]*mtt[15];
   mtt[16]=mtt[19] + mtt[17];
   mtt[16]=mtt[10]*mtt[16]*MMb;
   mtt[17]= - MMH*mtt[8];
   mtt[16]=mtt[17] + mtt[16];
   mtt[17]=mtt[4]*mtt[2];
   mtt[17]=3*mtt[17];
   mtt[18]=mtt[17]*mtt[7];
   mtt[18]=mtt[18] + mtt[23];
   mtt[18]=mtt[18]*MMb;
   mtt[17]=mtt[17]*mtt[5];
   mtt[19]=mtt[23] + mtt[8];
   mtt[17]= - mtt[17] + 1./2.*mtt[19];
   mtt[17]=MMt*mtt[17];
   mtt[19]=mtt[7] + mtt[14];
   mtt[19]=mtt[12] + mtt[5] + 1./2.*mtt[19];
   mtt[16]=1./8.*mtt[16] + mtt[17] - mtt[18] + 1./4.*mtt[19];
   mtt[17]=mtt[20] + mtt[21];
   mtt[16]=mtt[6]*mtt[17]*mtt[16];
   mtt[18]=MMZ*mtt[2];
   mtt[19]=mtt[14]*mtt[2];
   mtt[19]=mtt[18] - 1./4. + mtt[19];
   mtt[19]=1./4.*mtt[9] + 3*mtt[19] + mtt[23];
   mtt[19]=mtt[19]*mtt[20];
   mtt[23]= - 7./36.*mtt[9] - 1./36. + mtt[18];
   mtt[21]=mtt[23]*mtt[21];
   mtt[19]=mtt[19] + mtt[21];
   mtt[17]=mtt[2]*mtt[17];
   mtt[20]= - mtt[22] + 4./9. - 1./8.*mtt[20];
   mtt[20]=mtt[10]*mtt[20];
   mtt[17]=3./4.*mtt[17] + mtt[20];
   mtt[17]=mtt[13]*mtt[17];

      return  - 8./9. + 8./9.*mtt[9] + mtt[15] + mtt[16] + mtt[17] - 
      mtt[18] + 1./2.*mtt[19];
}
