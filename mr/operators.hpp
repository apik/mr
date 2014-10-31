//
// MR - 2-loop matching and 3-loop Running, including full 2-loop EW corrections
// Copyright (C) 2014 Andrey Pikelner <pikelner@theor.jinr.ru>
//
// This file is part of MR.
//
// MR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MR.  If not, see <http://www.gnu.org/licenses/>.
//

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
