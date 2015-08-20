#include <bb.hpp>
namespace mr
{
  long double bb<OS>::x01(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armbbbar[5], mbbbarret;

    armbbbar[1]=double(boson);
    armbbbar[2]=Tsil::A(MMb,mu2);
    armbbbar[3]=pow(MMb,-1);
    armbbbar[4]=armbbbar[2]*armbbbar[3];
    armbbbar[4]= - 1./3. + armbbbar[4];

    mbbbarret = 4*armbbbar[4]*armbbbar[1];
    return mbbbarret.real();
  }
} // namespace mr
