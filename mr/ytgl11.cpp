#include <tt.hpp>
std::complex<long double> tt::mygl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[13];

    mytgl[1]=pow(SW,-1);
    mytgl[2]=pow(MMH,-1);
    mytgl[3]=pow(MMW,-1);
    mytgl[4]=Tsil::B(MMH,MMt,MMt,mu2);
    mytgl[5]=Tsil::A(MMt,mu2);
    mytgl[6]=Tsil::A(MMH,mu2);
    mytgl[7]=pow(MMt,-1);
    mytgl[8]=std::real(Tsil::B(0,0,MMt,mu2));
   mytgl[9]= - 24*MMt - 28*mytgl[5];
   mytgl[9]=mytgl[2]*mytgl[9];
   mytgl[10]=97./4. + mytgl[8];
   mytgl[11]=pow(Pi,2);
   mytgl[9]= - 1./3.*mytgl[11] + 1./3.*mytgl[10] + 2*mytgl[4] + 
   mytgl[9];
   mytgl[9]=MMt*mytgl[9];
   mytgl[10]=3./2.*mytgl[6];
   mytgl[11]=1./4.*MMH + mytgl[10];
   mytgl[11]=mytgl[7]*mytgl[11];
   mytgl[12]= - 9*mytgl[2] - mytgl[7];
   mytgl[12]=mytgl[5]*mytgl[12];
   mytgl[11]=4*mytgl[12] + 6*mytgl[4] - 19./6. + 2*mytgl[8] + mytgl[11]
   ;
   mytgl[11]=mytgl[5]*mytgl[11];
   mytgl[12]= - 1./12. - mytgl[4];
   mytgl[12]=MMH*mytgl[12];
   mytgl[9]=mytgl[11] + mytgl[10] + mytgl[12] + mytgl[9];

      return mytgl[9]*mytgl[3]*pow(mytgl[1],2);
}
