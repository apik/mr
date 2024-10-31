#include <bb.hpp>
namespace mr
{
  double bb<MS>::x01(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> armbbos[5], mbbosret;

    armbbos[1]=double(boson);
    armbbos[2]=Tsil::A(mmb,mu2);
    armbbos[3]=pow(mmb,-1);
    armbbos[4]= - armbbos[3]*armbbos[2];
    armbbos[4]=1./3. + armbbos[4];

    mbbosret = 4*armbbos[4]*armbbos[1];
    return mbbosret.real();
  }
} // namespace mr
