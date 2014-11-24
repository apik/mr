#include <tt.hpp>
std::complex<long double>
tt::y01(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryutt[5], yuttret;

    aryutt[1]=double(boson);
    aryutt[2]=Tsil::A(MMt,mu2);
    aryutt[3]=pow(MMt,-1);
   aryutt[4]=aryutt[3]*aryutt[2];
   aryutt[4]= - 1./3. + aryutt[4];

      yuttret = 4*aryutt[4]*aryutt[1];
      return yuttret;
}
