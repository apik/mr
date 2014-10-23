#include <dr.hpp>
std::complex<long double> dr::dr10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mdr[18], mdrret;

    mdr[1]=pow(CW,-1);
    mdr[2]=pow(MMH,-1);
    mdr[3]=pow(MMZ,-1);
    mdr[4]=pow(SW,-1);
    mdr[5]=double(nH);
    mdr[6]=Tsil::A(MMt,mu2);
    mdr[7]=Tsil::A(MMH,mu2);
    mdr[8]=Tsil::A(MMZ,mu2);
    mdr[9]=Tsil::A(MMW,mu2);
    mdr[10]=1/( - MMW + MMH);
   mdr[11]=mdr[2]*MMZ;
   mdr[12]=mdr[11] - 1./8.;
   mdr[13]=mdr[8]*mdr[2];
   mdr[14]=mdr[9]*mdr[2];
   mdr[15]= - mdr[7] + mdr[9];
   mdr[15]=mdr[10]*mdr[15];
   mdr[14]=1./4.*mdr[15] + mdr[14] + 1./2.*mdr[13] + mdr[12];
   mdr[15]=mdr[8] - mdr[7];
   mdr[16]=mdr[6]*mdr[5];
   mdr[15]= - mdr[16] + mdr[9] + 1./2.*mdr[15];
   mdr[15]=1./4.*MMH + 3*mdr[15];
   mdr[17]=2*mdr[2];
   mdr[16]=mdr[17]*mdr[16];
   mdr[16]= - mdr[16] + 1./4.*mdr[5];
   mdr[17]=3*MMt;
   mdr[16]=mdr[16]*mdr[17];
   mdr[15]= - mdr[16] + 1./2.*mdr[15];
   mdr[16]=pow(mdr[4],2);
   mdr[17]=mdr[8] - mdr[9];
   mdr[17]=mdr[17]*mdr[16];
   mdr[17]=3./4.*mdr[17] - mdr[15];
   mdr[17]=mdr[3]*mdr[17];
   mdr[14]=3*mdr[14] + mdr[17];
   mdr[14]=mdr[14]*mdr[16];
   mdr[15]= - mdr[3]*mdr[15];
   mdr[12]=mdr[15] + 3./2.*mdr[13] + mdr[12];
   mdr[12]=mdr[12]*pow(mdr[1],2);

      mdrret =  - 2*mdr[11] + mdr[12] + mdr[14];
      return mdrret;
}
