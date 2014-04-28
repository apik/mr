// K. Ahnert and M. Mulansky, Odeint - Solving Ordinary Differential Equations in C++, AIP Conf. Proc. 1389, pp. 1586-1589 (2011);
// doi:http://dx.doi.org/10.1063/1.3637934
#ifndef __MR_HPP__
#define __MR_HPP__
// #include "tsil.hpp"
// #undef complex 
// #include <complex>
// #include <cmath>

#include "HH.hpp"
#include "WW.hpp"
#include "betaQCD.hpp"


template <typename T, typename U>
inline std::complex<T> operator*(std::complex<T> lhs, const U& rhs)
{
  return lhs *= rhs;
}

template <typename T, typename U>
inline std::complex<T> operator*(const U& rhs, std::complex<T> lhs)
{
  return lhs *= rhs;
}

template <typename T, typename U>
inline std::complex<T> operator/(std::complex<T> lhs, const U& rhs)
{
  return lhs /= rhs;
}


template <typename T, typename U>
inline std::complex<T> operator/(const U& lhs, std::complex<T> rhs)
{
  return std::complex<T>(lhs) /= rhs;
}

template <typename T, typename U>
inline std::complex<T> operator+(std::complex<T> lhs,const U& rhs)
{
  return lhs += rhs;
}

template <typename T, typename U>
inline std::complex<T> operator+(const U& rhs, std::complex<T> lhs)
{
  return lhs += rhs;
}

template <typename T, typename U>
inline std::complex<T> operator-(std::complex<T> lhs, const U& rhs)
{
  return lhs -= rhs;
}

template <typename T, typename U>
inline std::complex<T> operator-(const U& lhs, std::complex<T> rhs)
{
  return std::complex<T>(lhs) -= rhs;
}

#endif  // __MR_HPP__
