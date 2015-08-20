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

#include <tsil.hpp>

namespace mr
{
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

} // namespace mr
