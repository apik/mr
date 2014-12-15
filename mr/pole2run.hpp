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

#ifndef __POLE2RUN_HPP__
#define __POLE2RUN_HPP__

#include "sminput.hpp"
#include "betaSM.hpp"

class RunUpto
{
  
  long double al0;
  long double as0;
  long double mu0;

  // Constants at matching scale
  long double a1;
  long double a2;
  long double aS;
  long double ayt;
  long double alam;

  long double ms;

  CouplingsSM<3,3,3,3,-1,-1,3>* av;
public:
  RunUpto(OSinput oi, long double al_ = pdg2014::aMZ, long double as_ = pdg2014::asMZ, long double mu_ = pdg2014::MZ);
  
  long double lambda(long double mu);
  state_type operator()(long double mu);
};

#endif  // __POLE2RUN_HPP__
