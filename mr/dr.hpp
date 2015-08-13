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
  
    long double MMb,MMt, MMH, MMW, MMZ, mu2;
    long double SW, CW;
  
  
  public:
    dr()
    {
    }
  
    dr(long double,long double,long double,long double,long double);
  
    dr(OSinput, long double);
  
    void init();
  
    // 
    // Delta-r
    // 
    long double dr10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    long double dr11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    long double dr20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


    // Gaugeless limit
    long double drgl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    long double drgl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    long double drgl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  };


  template<>
  class dr<MS>
  {
  
    long double mmb,mmt, mmH, mmW, mmZ, mu2;
    long double s, c;
  
  
  public:
    dr()
    {
    }
  
    dr(long double,long double,long double,long double,long double);
  
    dr(MSinput, long double);
  
    void init();
  
    // 
    // Delta-r
    // 
    long double dr10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    long double dr11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    long double dr20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  };
} // namespace mr

#endif  //  __DR_HPP__
