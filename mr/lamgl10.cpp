#include <HH.hpp>
std::complex<long double> HH::lamgl10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mlamgl[10];

    mlamgl[1]=pow(SW,-1);
    mlamgl[2]=pow(MMW,-1);
    mlamgl[3]=Tsil::B(MMH,MMH,MMH,mu2);
    mlamgl[4]=Tsil::B(MMt,MMt,MMH,mu2);
    mlamgl[5]=pow(MMH,-1);
    mlamgl[6]=log(pow(mu2,-1)*MMt);
    mlamgl[7]=log(pow(mu2,-1)*MMH);
   mlamgl[8]=mlamgl[4] + 1./2. - mlamgl[6];
   mlamgl[9]=MMt*mlamgl[4]*mlamgl[5];
   mlamgl[8]=1./2.*mlamgl[8] - 2*mlamgl[9];
   mlamgl[8]=MMt*mlamgl[8];
   mlamgl[9]=9./8.*mlamgl[3] - 3./8.*mlamgl[7] + 7./8.;
   mlamgl[9]=MMH*mlamgl[9];
   mlamgl[8]=3*mlamgl[8] + mlamgl[9];

      return mlamgl[8]*mlamgl[2]*pow(mlamgl[1],2);
}
