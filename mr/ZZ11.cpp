#include <ZZ.hpp>
std::complex<long double> ZZ::m11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mZZ[16], mZZret;

    mZZ[1]=pow(CW,-1);
    mZZ[2]=pow(MMH,-1);
    mZZ[3]=pow(MMZ,-1);
    mZZ[4]=pow(SW,-1);
    mZZ[5]=Tsil::B(MMt,MMt,MMZ,mu2);
    mZZ[6]=Tsil::A(MMt,mu2);
    mZZ[7]=1/(4*MMt - MMZ);
   mZZ[8]=pow(mZZ[4],2);
   mZZ[9]=pow(mZZ[1],2);
   mZZ[10]=mZZ[8] + mZZ[9];
   mZZ[10]=mZZ[10]*mZZ[2];
   mZZ[11]=MMt*mZZ[10];
   mZZ[12]= - 7./9.*mZZ[9] + 64./9. + mZZ[8];
   mZZ[12]=mZZ[5]*mZZ[12];
   mZZ[13]=3*mZZ[8];
   mZZ[12]=mZZ[12] + 8*mZZ[11] - 43./9.*mZZ[9] + 64./9. - mZZ[13];
   mZZ[12]=MMt*mZZ[12];
   mZZ[11]= - mZZ[11] + 13./9.*mZZ[9] - 16./9. + mZZ[8];
   mZZ[14]=mZZ[13] - 64./3. + 25./3.*mZZ[9];
   mZZ[15]=mZZ[7]*mZZ[14];
   mZZ[10]= - 12*mZZ[10] + mZZ[15];
   mZZ[10]=mZZ[6]*mZZ[10];
   mZZ[13]=7./3.*mZZ[9] - 64./3. - mZZ[13];
   mZZ[13]=mZZ[5]*mZZ[13];
   mZZ[10]=4*mZZ[10] + 8*mZZ[11] + mZZ[13];
   mZZ[10]=mZZ[6]*mZZ[10];
   mZZ[10]=mZZ[10] + mZZ[12];
   mZZ[10]=mZZ[3]*mZZ[10];
   mZZ[11]= - 25./9.*mZZ[9] + 64./9. - mZZ[8];
   mZZ[11]=mZZ[7]*mZZ[11];
   mZZ[12]=mZZ[5]*mZZ[7];
   mZZ[13]=mZZ[14]*mZZ[12];
   mZZ[11]=4*mZZ[11] + mZZ[13];
   mZZ[11]=mZZ[6]*mZZ[11];
   mZZ[10]=mZZ[11] + mZZ[10];
   mZZ[8]=1./2.*mZZ[8] - 32./9. + 25./18.*mZZ[9];
   mZZ[9]=mZZ[7] - mZZ[12];
   mZZ[9]=MMZ*mZZ[9];
   mZZ[9]=mZZ[9] - mZZ[5] + 1;
   mZZ[8]=mZZ[8]*mZZ[9];

      mZZret = mZZ[8] + 2*mZZ[10];
      return mZZret;
}
