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
  std::complex<double> c2pp(TSIL_COMPLEX z)
  {
    return std::complex<double>(__real__ z, __imag__ z);
  }

  std::complex<double> csqrt(double z)
  {
    return sqrt(std::complex<double>(z,0));
  }

  std::complex<double> Li2(std::complex<double> z)
  {
    return TSIL_Dilog(*((TSIL_COMPLEX*) (&z)));
  }

  std::complex<double> Li3(std::complex<double> z)
  {
    return TSIL_Trilog(*((TSIL_COMPLEX*) (&z)));
  }

  std::complex<double> acc(std::complex<double> z)
  {    return z;  }

  std::complex<double> inv(std::complex<double> z)
  {    return std::complex<double>(1.,0)/z;  }

} // namespace mr
