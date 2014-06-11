#include <bb.hpp>
std::complex<long double> bb::mgl10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mbbgl[8];

    mbbgl[1]=pow(SW,-1);
    mbbgl[2]=pow(MMH,-1);
    mbbgl[3]=pow(MMW,-1);
    mbbgl[4]=log(pow(mu2,-1)*MMt);
    mbbgl[5]=log(pow(mu2,-1)*MMH);
   mbbgl[6]= - 5./2. + 3*mbbgl[4];
   mbbgl[7]=1 - mbbgl[4];
   mbbgl[7]=MMt*mbbgl[2]*mbbgl[7];
   mbbgl[6]=1./8.*mbbgl[6] + 3*mbbgl[7];
   mbbgl[6]=MMt*mbbgl[6];
   mbbgl[7]= - 1 + mbbgl[5];
   mbbgl[7]=MMH*mbbgl[7];
   mbbgl[6]=3./8.*mbbgl[7] + mbbgl[6];

      return mbbgl[6]*mbbgl[3]*pow(mbbgl[1],2);
}
