#include <bb.hpp>
std::complex<long double> bb::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbb[26];

    mbb[1]=pow(CW,-1);
    mbb[2]=pow(MMH,-1);
    mbb[3]=pow(SW,-1);
    mbb[4]=double(nH);
    mbb[5]=Tsil::A(MMt,mu2);
    mbb[6]=pow(MMZ,-1);
    mbb[7]=Tsil::A(MMb,mu2);
    mbb[8]=Tsil::B(MMH,MMb,MMb,mu2);
    mbb[9]=Tsil::B(MMZ,MMb,MMb,mu2);
    mbb[10]=pow(MMb,-1);
    mbb[11]=Tsil::B(MMW,MMt,MMb,mu2);
    mbb[12]=Tsil::A(MMH,mu2);
    mbb[13]=Tsil::A(MMZ,mu2);
    mbb[14]=Tsil::A(MMW,mu2);
   mbb[15]=mbb[13] - mbb[7];
   mbb[16]=MMt*mbb[11];
   mbb[17]=mbb[15] - mbb[16];
   mbb[18]= - 1./2.*mbb[9] - mbb[11];
   mbb[18]=MMZ*mbb[18];
   mbb[19]=mbb[14] - mbb[5];
   mbb[17]= - mbb[19] + mbb[18] - 1./2.*mbb[17];
   mbb[17]=mbb[10]*mbb[17];
   mbb[18]=mbb[14] + MMZ;
   mbb[20]=3*mbb[2];
   mbb[18]=mbb[20]*mbb[18];
   mbb[21]=mbb[9] - 3 + mbb[11];
   mbb[22]=mbb[13]*mbb[2];
   mbb[22]=3./2.*mbb[22];
   mbb[17]=1./2.*mbb[17] + mbb[22] + mbb[18] + 1./4.*mbb[21];
   mbb[18]=pow(mbb[3],2);
   mbb[17]=mbb[17]*mbb[18];
   mbb[21]=MMZ*mbb[2];
   mbb[22]=mbb[22] - 13./36. + mbb[21];
   mbb[23]=pow(mbb[1],2);
   mbb[22]=mbb[22]*mbb[23];
   mbb[17]=mbb[22] + mbb[17];
   mbb[22]=MMH*mbb[8];
   mbb[22]= - mbb[5] + mbb[22] - mbb[14];
   mbb[24]=1./4.*mbb[11];
   mbb[25]=mbb[24] + mbb[8];
   mbb[25]=mbb[25]*MMb;
   mbb[22]= - mbb[25] + 1./4.*mbb[22];
   mbb[20]=mbb[20]*mbb[4];
   mbb[25]=mbb[20]*mbb[5];
   mbb[25]=mbb[25] + mbb[24];
   mbb[25]=mbb[25]*MMt;
   mbb[16]=mbb[19] - mbb[16];
   mbb[19]=mbb[10]*MMt;
   mbb[19]=1./8.*mbb[19];
   mbb[16]=mbb[16]*mbb[19];
   mbb[19]=mbb[20]*MMb;
   mbb[19]=mbb[19] - 1./4.;
   mbb[19]=mbb[19]*mbb[7];
   mbb[16]=mbb[16] - 1./4.*mbb[12] + mbb[19] + mbb[25] + 1./2.*mbb[22];
   mbb[18]= - mbb[18] - mbb[23];
   mbb[16]=mbb[6]*mbb[16]*mbb[18];
   mbb[18]=1./9.*mbb[9];
   mbb[19]=1 - 5./8.*mbb[23];
   mbb[19]=mbb[18]*mbb[19];
   mbb[19]=mbb[19] + mbb[24];
   mbb[19]=MMZ*mbb[19];
   mbb[15]= - mbb[15]*mbb[23];
   mbb[15]=5./72.*mbb[15] + 1./9.*mbb[13] + 2./9.*mbb[7] + mbb[19];
   mbb[15]=mbb[10]*mbb[15];
   mbb[19]=2 + 17./8.*mbb[23];
   mbb[18]=mbb[19]*mbb[18];

      return  - 2./9. + mbb[15] + mbb[16] + 1./2.*mbb[17] + mbb[18] - 
      mbb[21];
}
