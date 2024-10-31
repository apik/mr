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

#ifndef __DR_HPP__
#define __DR_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"

namespace mr
{
  template<class T>
  class dr
  {
  };

  template<>
  class dr<OS>
  {
  
    double MMb,MMt, MMH, MMW, MMZ, mu2;
    double SW, CW;
  
  
  public:
    dr()
    {
    }
  
    dr(double,double,double,double,double);
  
    dr(OSinput, double);
  
    void init();
  
    // 
    // Delta-r
    // 
    double dr10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double dr11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double dr20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


    // Gaugeless limit
    double drgl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double drgl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double drgl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  };


  template<>
  class dr<MS>
  {
  
    double mmb,mmt, mmH, mmW, mmZ, mu2;
    double s, c;
  
  
  public:
    dr()
    {
    }
  
    dr(double,double,double,double,double);
  
    dr(MSinput, double);
  
    void init();
  
    // 
    // Delta-r
    // 
    double dr10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double dr11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double dr20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  };
} // namespace mr

#endif  //  __DR_HPP__
