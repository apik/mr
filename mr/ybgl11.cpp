#include <bb.hpp>
std::complex<long double> bb::mygl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mybgl[11];

    mybgl[1]=pow(SW,-1);
    mybgl[2]=pow(MMW,-1);
    mybgl[3]=log(pow(mu2,-1)*MMb);
    mybgl[4]=log(pow(mu2,-1)*MMt);
    mybgl[5]=log(MMt);
    mybgl[6]=log(MMb);
   mybgl[7]=1./4.*mybgl[3];
   mybgl[8]=pow(Pi,2);
   mybgl[9]=mybgl[6] + mybgl[5];
   mybgl[9]=5 - 3./2.*mybgl[9];
   mybgl[9]=mybgl[4]*mybgl[9];
   mybgl[9]=1./3.*mybgl[8] + mybgl[7] + 11./12. + mybgl[9];
   mybgl[9]=MMt*mybgl[9];
   mybgl[10]=2 - mybgl[6];
   mybgl[10]=mybgl[3]*mybgl[10];
   mybgl[8]= - 1./12.*mybgl[8] - 3./2. + mybgl[10];
   mybgl[8]=MMb*mybgl[8];
   mybgl[7]= - 1./3. + mybgl[7];
   mybgl[7]=MMH*mybgl[7];
   mybgl[7]=mybgl[7] + 5*mybgl[8] + mybgl[9];

      return mybgl[7]*mybgl[2]*pow(mybgl[1],2);
}
