#include <tt.hpp>
std::complex<long double> tt::drgl10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mdrgl[8];

    mdrgl[1]=pow(SW,-1);
    mdrgl[2]=pow(MMW,-1);
    mdrgl[3]=Tsil::A(MMH,mu2);
    mdrgl[4]=Tsil::A(MMt,mu2);
    mdrgl[5]=pow(MMH,-1);
   mdrgl[6]=MMt + mdrgl[3];
   mdrgl[7]=mdrgl[5]*MMt;
   mdrgl[7]=1./2. - 2*mdrgl[7];
   mdrgl[7]=mdrgl[4]*mdrgl[7];
   mdrgl[6]=1./4.*mdrgl[6] + mdrgl[7];
   mdrgl[6]=3*mdrgl[6] - 1./8.*MMH;

      return mdrgl[6]*mdrgl[2]*pow(mdrgl[1],2);
}
