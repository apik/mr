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

#ifndef __MFITTER_HPP__
#define __MFITTER_HPP__ 

#include <vector>
#include "sminput.hpp"
#include "Minuit/FCNBase.h"


std::vector<long double> observables(long double,
                                     long double,
                                     long double,
                                     long double,
                                     long double,
                                     long double,
                                     long double);


class MfitterFcn : public FCNBase {
  
public:
  
  MfitterFcn(const std::vector<double>& meas,
             // const std::vector<double>& pos,
             const std::vector<double>& mvar,
             double scale) : theMeasurements(meas),
                             // thePositions(pos),
                             theMVariances(mvar),
                             theScale(scale),
                             theErrorDef(1.) {}

  ~MfitterFcn() {}

  virtual double up() const {return theErrorDef;}
  virtual double operator()(const std::vector<double>&) const;
  
  std::vector<double> measurements() const {return theMeasurements;}
  // std::vector<double> positions() const {return thePositions;}
  std::vector<double> variances() const {return theMVariances;}

  void setErrorDef(double def) {theErrorDef = def;}

private:

  
  std::vector<double> theMeasurements;
  // std::vector<double> thePositions;
  std::vector<double> theMVariances;
  double theErrorDef;
  double theScale;
};


#endif  // __MFITTER_HPP__ 
