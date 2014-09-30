#ifndef __OPERATORS_HPP__
#define __OPERATORS_HPP__

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

template <typename T>
inline std::complex<T> Re(std::complex<T> z)
{
  return std::complex<T>(z.real(),0);
}

#endif  // __OPERATORS_HPP__
