#include <HH.hpp>
std::complex<long double> HH::mgl10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHHgl[10];

    mHHgl[1]=pow(SW,-1);
    mHHgl[2]=pow(MMH,-1);
    mHHgl[3]=pow(MMW,-1);
    mHHgl[4]=Tsil::B(MMH,MMH,MMH,mu2);
    mHHgl[5]=Tsil::B(MMt,MMt,MMH,mu2);
    mHHgl[6]=log(pow(mu2,-1)*MMt);
    mHHgl[7]=log(pow(mu2,-1)*MMH);
   mHHgl[8]=1./8.*mHHgl[7] + 3./8.*mHHgl[4];
   mHHgl[8]=MMH*mHHgl[8];
   mHHgl[9]= - mHHgl[6] + 1 - mHHgl[5];
   mHHgl[9]=MMt*mHHgl[2]*mHHgl[9];
   mHHgl[9]=1./2.*mHHgl[5] + 2*mHHgl[9];
   mHHgl[9]=MMt*mHHgl[9];
   mHHgl[8]=mHHgl[9] + mHHgl[8];

      return 3*mHHgl[8]*mHHgl[3]*pow(mHHgl[1],2);
}
