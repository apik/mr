#include <tt.hpp>
namespace mr
{
  double tt<OS>::x01(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> armttbar[5], mttbarret;

    armttbar[1]=double(boson);
    armttbar[2]=Tsil::A(MMt,mu2);
    armttbar[3]=pow(MMt,-1);
    armttbar[4]=armttbar[2]*armttbar[3];
    armttbar[4]= - 1./3. + armttbar[4];

    mttbarret = 4*armttbar[4]*armttbar[1];
    return mttbarret.real();
  }
} // namespace mr
