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

  Tsil* protWHHWW;
  Tsil* protWHZWW;
  Tsil* protWZZWW;
  Tsil* protWWHHH;
  Tsil* protWWHZZ;
  Tsil* protWWZZH;
  Tsil* protWtZ00;
  Tsil* protW0HWW;
  Tsil* protW0Htt;
  Tsil* protW0ZWW;
  Tsil* protW0Ztt;
  Tsil* protW0Z00;
  Tsil* prot0WW0W;
  Tsil* prot0Wt0t;
  Tsil* prot0W0Z0;
  Tsil* prot00Wt0;
  Tsil* prot00W00;
  Tsil* prot00ttZ;
  Tsil* prot00tt0;
  Tsil* prot0000Z;
  Tsil* prot00000;
  // ZZ part

  Tsil* protZHHZZ;
  Tsil* protZZHHH;
  Tsil* protZWHWW;
  Tsil* prottZtHt;
  Tsil* protWWWWH;
  Tsil* protWWWWZ;
  Tsil* protWWWW0;
  Tsil* protWtWt0;
  Tsil* protW0W0t;
  Tsil* protW0W00;
  Tsil* protttttH;
  Tsil* protttttZ;
  Tsil* prottttt0;
  Tsil* prott0t0W;
  // Tsil* prot0000Z;
  Tsil* prot0000W;
  // Tsil* prot00000;
  Tsil* protWZWHW;
  
  TsilSTU* protHZ00;


  TsilSTU* protWH0H;
  TsilSTU* protWZ0Z;
  TsilSTU* protHW00;
  
  TsilST* protos[24];
public:
  alphaGF()
  {
  }

  alphaGF(long double,long double,long double,long double,long double);

  alphaGF(OSinput, long double);

  void init();
  
 
  std::complex<long double> a10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  std::complex<long double> a11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

  std::complex<long double> a20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
};


#endif  // __ALPHAGF_HPP__
