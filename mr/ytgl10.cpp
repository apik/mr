#include <tt.hpp>
std::complex<long double> tt::mygl10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[9];

    mytgl[1]=pow(SW,-1);
    mytgl[2]=pow(MMW,-1);
    mytgl[3]=Tsil::B(MMH,MMt,MMt,mu2);
    mytgl[4]=Tsil::B(0,0,MMt,mu2);
    mytgl[5]=log(pow(mu2,-1)*MMt);
    mytgl[6]=log(pow(mu2,-1)*MMH);
   mytgl[7]= - mytgl[3] + 3./2. - mytgl[6];
   mytgl[7]=MMH*mytgl[7];
   mytgl[8]= - mytgl[5] + 1./4.*mytgl[4] + 1./4. + mytgl[3];
   mytgl[8]=MMt*mytgl[8];
   mytgl[7]=1./4.*mytgl[7] + mytgl[8];

      return 1./2.*mytgl[7]*mytgl[2]*pow(mytgl[1],2);
}
