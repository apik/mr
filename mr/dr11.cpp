#include <dr.hpp>
std::complex<long double> dr::dr11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mdr[12], mdrret;

    mdr[1]=pow(CW,-1);
    mdr[2]=pow(MMH,-1);
    mdr[3]=pow(MMZ,-1);
    mdr[4]=pow(SW,-1);
    mdr[5]=Tsil::I2(0,0,MMt,mu2);
    mdr[6]=Tsil::A(MMt,mu2);
    mdr[7]=pow(MMt,-1);
    mdr[8]=Tsil::Aeps(MMt,mu2);
   mdr[9]= - mdr[7] + 12*mdr[2];
   mdr[10]=2*mdr[6];
   mdr[9]=mdr[9]*mdr[10];
   mdr[9]=mdr[9] - 5;
   mdr[9]=mdr[9]*mdr[10];
   mdr[10]=mdr[5] - mdr[8];
   mdr[11]=MMt*mdr[2];
   mdr[11]= - 37./2. + 64*mdr[11];
   mdr[11]=mdr[11]*MMt;
   mdr[9]= - mdr[9] + mdr[11] - 4*mdr[10];
   mdr[10]=pow(mdr[1],2);
   mdr[11]=pow(mdr[4],2);
   mdr[10]=mdr[10] + mdr[11];

      mdrret = mdr[10]*mdr[9]*mdr[3];
      return mdrret;
}
