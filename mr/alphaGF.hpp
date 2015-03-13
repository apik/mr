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

class alphaGF
{

  long double MMb, MMt, MMH, MMW, MMZ, mu2;
  long double SW, CW;

  Tsil* WprotWHHWW;
  Tsil* WprotWHZWW;
  Tsil* WprotWZZWW;
  Tsil* WprotWWHHH;
  Tsil* WprotWWHZZ;
  Tsil* WprotWWZZH;
  Tsil* WprotWtZ00;
  Tsil* WprotW0HWW;
  Tsil* WprotW0Htt;
  Tsil* WprotW0ZWW;
  Tsil* WprotW0Ztt;
  Tsil* WprotW0Z00;
  Tsil* Wprot0WW0W;
  Tsil* Wprot0Wt0t;
  Tsil* Wprot0W0Z0;
  Tsil* Wprot00Wt0;
  Tsil* Wprot00W00;
  Tsil* Wprot00ttZ;
  Tsil* Wprot00tt0;
  Tsil* Wprot0000Z;
  Tsil* Wprot00000;
  TsilSTU* WprotHW00;
  TsilSTU* WprotWH0H;
  TsilSTU* WprotWZ0Z;

  
  // ZZ part
  Tsil* ZprotZHHZZ;
  Tsil* ZprotZZHHH;
  Tsil* ZprotZWHWW;
  Tsil* ZprottZtHt;
  Tsil* ZprotWWWWH;
  Tsil* ZprotWWWWZ;
  Tsil* ZprotWWWW0;
  Tsil* ZprotWtWt0;
  Tsil* ZprotW0W0t;
  Tsil* ZprotW0W00;
  Tsil* ZprotttttH;
  Tsil* ZprotttttZ;
  Tsil* Zprottttt0;
  Tsil* Zprott0t0W;
  Tsil* Zprot0000Z;
  Tsil* Zprot0000W;
  Tsil* Zprot00000;
  Tsil* ZprotWZWHW;
  TsilSTU* ZprotHZ00;



  
  
  TsilSTU* Wprotos[24];
  TsilSTU* Zprotos[19];
  
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


#endif  // __ALPHAGF_HPP__
