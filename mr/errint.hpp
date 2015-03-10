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

#ifndef __ERRINT_HPP__
#define __ERRINT_HPP__

#include "boost/numeric/interval.hpp"
#include "sminput.hpp"
// interval type for Long Double numbers
typedef boost::numeric::interval<long double> LDI;

typedef OSinputTemplate<LDI> OSinputErr;

template <class T>
boost::numeric::interval<T> we(T mean, T err)
{
  return boost::numeric::interval<T>(mean - err, mean + err);
}

class WWerr
{
public:
  WWerr(OSinputErr, long double);
};

#endif  // __ERRINT_HPP__
