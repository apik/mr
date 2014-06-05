#include <tsil.hpp>

std::complex<long double> c2pp(TSIL_COMPLEX z)
{
  return std::complex<long double>(__real__ z, __imag__ z);
}

std::complex<long double> csqrt(long double z)
{
  return sqrt(std::complex<long double>(z,0));
}

std::complex<long double> Li2(std::complex<long double> z)
{
  return TSIL_Dilog(z.real()+1.0I*z.imag());
}

std::complex<long double> Li3(std::complex<long double> z)
{
  return TSIL_Trilog(z.real()+1.0I*z.imag());
}

std::complex<long double> acc(std::complex<long double> z)
{    return z;  }  

std::complex<long double> inv(std::complex<long double> z)
{    return std::complex<long double>(1.,0)/z;  }  
