#include <tt.hpp>
std::complex<long double> tt::mgl10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mttgl[10];

    mttgl[1]=Tsil::B(MMH,MMt,MMt,mu2);
    mttgl[2]=pow(SW,-1);
    mttgl[3]=pow(MMW,-1);
    mttgl[4]=Tsil::A(MMH,mu2);
    mttgl[5]=Tsil::A(MMt,mu2);
    mttgl[6]=pow(MMH,-1);
    mttgl[7]=std::real(Tsil::B(0,0,MMt,mu2));
   mttgl[8]=1./2.*mttgl[1];
   mttgl[9]=mttgl[6]*mttgl[5];
   mttgl[9]=mttgl[8] + 1./8.*mttgl[7] - 3*mttgl[9];
   mttgl[9]=MMt*mttgl[9];
   mttgl[8]= - MMH*mttgl[8];
   mttgl[8]=mttgl[8] + mttgl[4] + mttgl[5];
   mttgl[8]=1./4.*mttgl[8] + mttgl[9];

      return mttgl[8]*mttgl[3]*pow(mttgl[2],2);
}
