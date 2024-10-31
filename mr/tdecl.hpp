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

#ifndef __TDECL_HPP__
#define __TDECL_HPP__


// Type for couplings 
typedef double Rt;

// Type for masses
typedef double MRt;

// Couplings set
typedef std::vector<Rt> SMCouplings;


//
// Tolerance class for root-finding algorithms
//
class Tolerance {
public:
  Tolerance(double eps) :
    _eps(eps) {
  }
  bool operator()(double a, double b) {
    return (std::abs(b - a) <= _eps);
  }
private:
  double _eps;
};


#endif  // __TDECL_HPP__
