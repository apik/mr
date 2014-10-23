#include <dr.hpp>
std::complex<long double> dr::dr11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mdr[13], mdrret;

    mdr[1]=double(nH);
    mdr[2]=pow(CW,-1);
    mdr[3]=pow(MMH,-1);
    mdr[4]=pow(MMZ,-1);
    mdr[5]=pow(SW,-1);
    mdr[6]=Tsil::I2(0,0,MMt,mu2);
    mdr[7]=Tsil::A(MMt,mu2);
    mdr[8]=pow(MMt,-1);
    mdr[9]=Tsil::Aeps(MMt,mu2);
   mdr[10]= - mdr[8] + 12*mdr[3];
   mdr[11]=2*mdr[7];
   mdr[10]=mdr[10]*mdr[11];
   mdr[10]=mdr[10] - 5;
   mdr[10]=mdr[10]*mdr[11];
   mdr[11]=mdr[9] - mdr[6];
   mdr[12]=MMt*mdr[3];
   mdr[12]= - 37./2. + 64*mdr[12];
   mdr[12]=mdr[12]*MMt;
   mdr[10]= - mdr[10] + mdr[12] + 4*mdr[11];
   mdr[11]=pow(mdr[2],2);
   mdr[12]=pow(mdr[5],2);
   mdr[11]=mdr[11] + mdr[12];

      mdrret = mdr[11]*mdr[10]*mdr[4]*mdr[1];
      return mdrret;
}
