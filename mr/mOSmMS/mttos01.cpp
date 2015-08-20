#include <tt.hpp>
namespace mr
{
  long double tt<MS>::x01(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armttos[5], mttosret;

    armttos[1]=double(boson);
    armttos[2]=Tsil::A(mmt,mu2);
    armttos[3]=pow(mmt,-1);
    armttos[4]= - armttos[3]*armttos[2];
    armttos[4]=1./3. + armttos[4];

    mttosret = 4*armttos[4]*armttos[1];
    return mttosret.real();
  }
} // namespace mr
