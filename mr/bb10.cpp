#include <bb.hpp>
std::complex<long double> bb::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbb[19];

    mbb[1]=pow(CW,-1);
    mbb[2]=pow(MMH,-1);
    mbb[3]=pow(MMZ,-1);
    mbb[4]=pow(SW,-1);
    mbb[5]=Tsil::A(MMH,mu2);
    mbb[6]=Tsil::A(MMZ,mu2);
    mbb[7]=Tsil::A(MMW,mu2);
    mbb[8]=Tsil::A(MMt,mu2);
    mbb[9]=Tsil::A(MMb,mu2);
    mbb[10]=pow(MMb,-1);
    mbb[11]=1/(MMt - MMW);
   mbb[12]=1./2.*mbb[6] + mbb[7];
   mbb[12]=mbb[2]*mbb[12];
   mbb[13]=MMZ*mbb[2];
   mbb[12]=mbb[13] - 1./8. + mbb[12];
   mbb[14]=mbb[5] + mbb[8];
   mbb[15]=3*mbb[2];
   mbb[16]=mbb[15]*mbb[8];
   mbb[16]=mbb[16] - 1./16.;
   mbb[16]=mbb[16]*MMt;
   mbb[14]=mbb[16] - 3./8.*mbb[14];
   mbb[16]= - mbb[3]*mbb[14];
   mbb[17]=3./8.*mbb[11];
   mbb[18]=mbb[7] - mbb[8];
   mbb[18]=mbb[18]*mbb[11];
   mbb[18]=mbb[18] + 1;
   mbb[18]=MMZ*mbb[18];
   mbb[19]=mbb[7] + mbb[18];
   mbb[19]=mbb[19]*mbb[17];
   mbb[12]= - 1./6. + mbb[19] + 3./2.*mbb[12] + mbb[16];
   mbb[12]=mbb[12]*pow(mbb[4],2);
   mbb[15]=mbb[6]*mbb[15];
   mbb[15]= - 31./36. + mbb[15];
   mbb[15]=1./2.*mbb[15] + mbb[13];
   mbb[14]= - 1./6.*mbb[6] - mbb[14];
   mbb[14]=mbb[3]*mbb[14];
   mbb[14]=1./2.*mbb[15] + mbb[14];
   mbb[14]=mbb[14]*pow(mbb[1],2);
   mbb[15]=mbb[3]*mbb[6];
   mbb[16]=mbb[10]*mbb[9];
   mbb[15]=mbb[15] - mbb[16];
   mbb[16]= - mbb[18]*mbb[17];

      return mbb[12] - mbb[13] + mbb[14] - 1./3.*mbb[15] + mbb[16];
}
