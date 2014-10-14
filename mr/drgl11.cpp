#include <dr.hpp>
std::complex<long double> dr::drgl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mdrgl[8], mdrglret;

    mdrgl[1]=pow(SW,-1);
    mdrgl[2]=pow(MMH,-1);
    mdrgl[3]=pow(MMW,-1);
    mdrgl[4]=Tsil::A(MMt,mu2);
    mdrgl[5]=pow(MMt,-1);
   mdrgl[6]=pow(Pi,2);
   mdrgl[7]=MMt*mdrgl[2];
   mdrgl[6]=64*mdrgl[7] - 29./2. + 2./3.*mdrgl[6];
   mdrgl[6]=MMt*mdrgl[6];
   mdrgl[7]= - 8*mdrgl[2] + mdrgl[5];
   mdrgl[7]=mdrgl[4]*mdrgl[7];
   mdrgl[7]=1 + mdrgl[7];
   mdrgl[7]=mdrgl[4]*mdrgl[7];
   mdrgl[6]=mdrgl[6] + 6*mdrgl[7];

      mdrglret = mdrgl[6]*mdrgl[3]*pow(mdrgl[1],2);
      return mdrglret;
}
