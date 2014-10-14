#include <dr.hpp>
std::complex<long double> dr::dr10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mdr[16], mdrret;

    mdr[1]=pow(CW,-1);
    mdr[2]=pow(MMH,-1);
    mdr[3]=pow(MMZ,-1);
    mdr[4]=pow(SW,-1);
    mdr[5]=Tsil::A(MMH,mu2);
    mdr[6]=Tsil::A(MMZ,mu2);
    mdr[7]=Tsil::A(MMW,mu2);
    mdr[8]=Tsil::A(MMt,mu2);
    mdr[9]=1/( - MMW + MMH);
   mdr[10]=2*mdr[2];
   mdr[10]=mdr[10]*mdr[8];
   mdr[10]=mdr[10] - 1./4.;
   mdr[10]=mdr[10]*MMt;
   mdr[11]=mdr[7] - mdr[8];
   mdr[10]=mdr[10] + 1./2.*mdr[11];
   mdr[10]= - 3./4.*mdr[5] + 1./8.*MMH + 3*mdr[10];
   mdr[10]=mdr[10]*mdr[3];
   mdr[11]= - mdr[2] + 1./2.*mdr[3];
   mdr[12]=3./2.*mdr[6];
   mdr[11]=mdr[11]*mdr[12];
   mdr[10]=mdr[10] + mdr[11];
   mdr[11]= - mdr[7] + mdr[6];
   mdr[12]=pow(mdr[4],2);
   mdr[11]=mdr[12]*mdr[3]*mdr[11];
   mdr[13]=mdr[7] - mdr[5];
   mdr[13]=mdr[9]*mdr[13];
   mdr[11]=mdr[11] + mdr[13];
   mdr[13]=mdr[2]*MMZ;
   mdr[14]=mdr[13] - 1./8.;
   mdr[15]=mdr[7]*mdr[2];
   mdr[15]=mdr[15] + mdr[14];
   mdr[11]=3*mdr[15] - mdr[10] + 3./4.*mdr[11];
   mdr[11]=mdr[12]*mdr[11];
   mdr[10]=mdr[14] - mdr[10];
   mdr[10]=mdr[10]*pow(mdr[1],2);

      mdrret = mdr[10] + mdr[11] - 2*mdr[13];
      return mdrret;
}
