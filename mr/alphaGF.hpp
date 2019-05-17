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

#ifndef __ALPHAGF_HPP__
#define __ALPHAGF_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"
#include "base.hpp"

namespace mr
{
  class alphaGF
  {

    long double MMb, MMt, MMH, MMW, MMZ, mu2;
    long double SW, CW;

    std::unique_ptr<Tsil> WprotWHHWW;
    std::unique_ptr<Tsil> WprotWHZWW;
    std::unique_ptr<Tsil> WprotWZZWW;
    std::unique_ptr<Tsil> WprotWWHHH;
    std::unique_ptr<Tsil> WprotWWHZZ;
    std::unique_ptr<Tsil> WprotWWZZH;
    std::unique_ptr<Tsil> WprotWtZ00;
    std::unique_ptr<Tsil> WprotW0HWW;
    std::unique_ptr<Tsil> WprotW0Htt;
    std::unique_ptr<Tsil> WprotW0ZWW;
    std::unique_ptr<Tsil> WprotW0Ztt;
    std::unique_ptr<Tsil> WprotW0Z00;
    std::unique_ptr<Tsil> Wprot0WW0W;
    std::unique_ptr<Tsil> Wprot0Wt0t;
    std::unique_ptr<Tsil> Wprot0W0Z0;
    std::unique_ptr<Tsil> Wprot00Wt0;
    std::unique_ptr<Tsil> Wprot00W00;
    std::unique_ptr<Tsil> Wprot00ttZ;
    std::unique_ptr<Tsil> Wprot00tt0;
    std::unique_ptr<Tsil> Wprot0000Z;
    std::unique_ptr<Tsil> Wprot00000;
    std::unique_ptr<TsilSTU> WprotHW00;
    std::unique_ptr<TsilSTU> WprotWH0H;
    std::unique_ptr<TsilSTU> WprotWZ0Z;

  
    // ZZ part
    std::unique_ptr<Tsil> ZprotZHHZZ;
    std::unique_ptr<Tsil> ZprotZZHHH;
    std::unique_ptr<Tsil> ZprotZWHWW;
    std::unique_ptr<Tsil> ZprottZtHt;
    std::unique_ptr<Tsil> ZprotWWWWH;
    std::unique_ptr<Tsil> ZprotWWWWZ;
    std::unique_ptr<Tsil> ZprotWWWW0;
    std::unique_ptr<Tsil> ZprotWtWt0;
    std::unique_ptr<Tsil> ZprotW0W0t;
    std::unique_ptr<Tsil> ZprotW0W00;
    std::unique_ptr<Tsil> ZprotttttH;
    std::unique_ptr<Tsil> ZprotttttZ;
    std::unique_ptr<Tsil> Zprottttt0;
    std::unique_ptr<Tsil> Zprott0t0W;
    std::unique_ptr<Tsil> Zprot0000Z;
    std::unique_ptr<Tsil> Zprot0000W;
    std::unique_ptr<Tsil> Zprot00000;
    std::unique_ptr<Tsil> ZprotWZWHW;
    std::unique_ptr<TsilSTU> ZprotHZ00;
  
  public:
    alphaGF()
    {
    }

    alphaGF(long double,long double,long double,long double,long double);

    alphaGF(OSinput, long double);

    void init();
  
 
    long double a10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
    long double a11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

    long double a20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  };
} // namespace mr

#endif  // __ALPHAGF_HPP__
