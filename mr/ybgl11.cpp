#include <bb.hpp>
std::complex<long double> bb::mygl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mybgl[12];

    mybgl[1]=pow(SW,-1);
    mybgl[2]=pow(MMW,-1);
    mybgl[3]=Tsil::I2(0,0,MMt,mu2);
    mybgl[4]=Tsil::Aeps(MMt,mu2);
    mybgl[5]=Tsil::Aeps(MMb,mu2);
    mybgl[6]=log(pow(mu2,-1)*MMb);
    mybgl[7]=log(pow(mu2,-1)*MMt);
   mybgl[8]=1./4.*mybgl[6];
   mybgl[9]=pow(Pi,2);
   mybgl[9]= - 109./4. - mybgl[9];
   mybgl[10]= - 7./2.*mybgl[7] + 13 - 3./2.*mybgl[6];
   mybgl[10]=mybgl[7]*mybgl[10];
   mybgl[9]=mybgl[10] + 1./3.*mybgl[9] + mybgl[8];
   mybgl[9]=MMt*mybgl[9];
   mybgl[10]=1./2.*MMb;
   mybgl[11]= - mybgl[6]*mybgl[10];
   mybgl[11]=MMb + mybgl[11];
   mybgl[11]=mybgl[6]*mybgl[11];
   mybgl[10]=mybgl[5] - mybgl[10] + mybgl[11];
   mybgl[11]=mybgl[3] - mybgl[4];
   mybgl[8]= - 1./3. + mybgl[8];
   mybgl[8]=MMH*mybgl[8];
   mybgl[8]=mybgl[8] + mybgl[9] + 5*mybgl[10] - 4*mybgl[11];

      return mybgl[8]*mybgl[2]*pow(mybgl[1],2);
}
