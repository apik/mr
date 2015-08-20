#include <bb.hpp>
namespace mr
{
  long double bb<OS>::y01(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> aryubb[5], yubbret;

    aryubb[1]=double(boson);
    aryubb[2]=Tsil::A(MMb,mu2);
    aryubb[3]=pow(MMb,-1);
    aryubb[4]=aryubb[2]*aryubb[3];
    aryubb[4]= - 1./3. + aryubb[4];

    yubbret = 4*aryubb[4]*aryubb[1];
    return yubbret.real();
  }
} // namespace mr
