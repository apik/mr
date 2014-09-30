#include <bb.hpp>
std::complex<long double> bb::mygl10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mybgl[6];

    mybgl[1]=pow(SW,-1);
    mybgl[2]=pow(MMW,-1);
    mybgl[3]=log(pow(mu2,-1)*MMt);
   mybgl[4]=MMH + MMt;
   mybgl[5]=mybgl[3]*MMt;
   mybgl[4]=1./2.*mybgl[4] - 3*mybgl[5];

      return 1./8.*mybgl[4]*mybgl[2]*pow(mybgl[1],2);
}
