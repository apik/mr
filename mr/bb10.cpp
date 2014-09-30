#include <bb.hpp>
std::complex<long double> bb::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbb[19];

    mbb[1]=pow(CW,-1);
    mbb[2]=pow(MMZ,-1);
    mbb[3]=pow(MMt,-1);
    mbb[4]=pow(SW,-1);
    mbb[5]=Tsil::A(MMZ,mu2);
    mbb[6]=Tsil::A(MMW,mu2);
    mbb[7]=Tsil::A(MMt,mu2);
    mbb[8]=Tsil::A(MMb,mu2);
    mbb[9]=pow(MMb,-1);
    mbb[10]=1/(MMt - MMW);
   mbb[11]=pow(mbb[1],2);
   mbb[12]=1./2.*mbb[11];
   mbb[13]=pow(mbb[4],2);
   mbb[14]=mbb[12] - 1 + 3./2.*mbb[13];
   mbb[14]=MMZ*mbb[14];
   mbb[15]=mbb[11] + mbb[13];
   mbb[16]=mbb[5]*mbb[15];
   mbb[17]=mbb[13]*mbb[6];
   mbb[14]=3./4.*mbb[16] + 3./2.*mbb[17] + mbb[14];
   mbb[14]=mbb[3]*mbb[14];
   mbb[16]=mbb[6] - mbb[7];
   mbb[16]=mbb[10]*mbb[16];
   mbb[16]=mbb[16] + 1;
   mbb[18]=mbb[13] - 1;
   mbb[16]=MMZ*mbb[18]*mbb[16];
   mbb[16]=mbb[16] + mbb[17];
   mbb[17]=3./8.*mbb[10];
   mbb[16]=mbb[17]*mbb[16];
   mbb[17]=1./16.*MMt - 9./4.*mbb[7];
   mbb[15]=mbb[15]*mbb[17];
   mbb[12]= - 1 - mbb[12];
   mbb[12]=mbb[5]*mbb[12];
   mbb[12]= - 1./6. + 1./3.*mbb[12] + mbb[15];
   mbb[12]=mbb[2]*mbb[12];
   mbb[15]=mbb[9]*mbb[8];

      return  - 31./144.*mbb[11] + mbb[12] - 3./16.*mbb[13] + mbb[14]
       + 1./3.*mbb[15] + mbb[16];
}
