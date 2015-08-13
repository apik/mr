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

#include <alphaGF.hpp>
#include "timer.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif


alphaGF::alphaGF(OSinput sm, long double mu2_)
{
  MMb = sm.MMb();
  MMW = sm.MMW();
  MMZ = sm.MMZ();
  MMH = sm.MMH();
  MMt = sm.MMt();
  mu2 = mu2_;

  init();
}


void alphaGF::init()
{

  CW = sqrt(MMW/MMZ);
  SW = sqrt(1-MMW/MMZ);

  Wprotos[0] =  WprotWHHWW = new Tsil(MMW, MMH, MMH, MMW, MMW, mu2);
  Wprotos[1] =  WprotWHZWW = new Tsil(MMW, MMH, MMZ, MMW, MMW, mu2);
  Wprotos[2] =  WprotWZZWW = new Tsil(MMW, MMZ, MMZ, MMW, MMW, mu2);
  Wprotos[3] =  WprotWWHHH = new Tsil(MMW, MMW, MMH, MMH, MMH, mu2);
  Wprotos[4] =  WprotWWHZZ = new Tsil(MMW, MMW, MMH, MMZ, MMZ, mu2);
  Wprotos[5] =  WprotWWZZH = new Tsil(MMW, MMW, MMZ, MMZ, MMH, mu2);
  Wprotos[6] =  WprotWtZ00 = new Tsil(MMW, MMt, MMZ,   0,   0, mu2);
  Wprotos[7] =  WprotW0HWW = new Tsil(MMW,   0, MMH, MMW, MMW, mu2);
  Wprotos[8] =  WprotW0Htt = new Tsil(MMW,   0, MMH, MMt, MMt, mu2);
  Wprotos[9] =  WprotW0ZWW = new Tsil(MMW,   0, MMZ, MMW, MMW, mu2);
  Wprotos[10] =  WprotW0Ztt = new Tsil(MMW,   0, MMZ, MMt, MMt, mu2);
  Wprotos[11] =  WprotW0Z00 = new Tsil(MMW,   0, MMZ,   0,   0, mu2);
  Wprotos[12] =  Wprot0WW0W = new Tsil(  0, MMW, MMW,   0, MMW, mu2);
  Wprotos[13] =  Wprot0Wt0t = new Tsil(  0, MMW, MMt,   0, MMt, mu2);
  Wprotos[14] =  Wprot0W0Z0 = new Tsil(  0, MMW,   0, MMZ,   0, mu2);
  Wprotos[15] =  Wprot00Wt0 = new Tsil(  0,   0, MMW, MMt,   0, mu2);
  Wprotos[16] =  Wprot00W00 = new Tsil(  0,   0, MMW,   0,   0, mu2);
  Wprotos[17] =  Wprot00ttZ = new Tsil(  0,   0, MMt, MMt, MMZ, mu2);
  Wprotos[18] =  Wprot00tt0 = new Tsil(  0,   0, MMt, MMt,   0, mu2);
  Wprotos[19] =  Wprot0000Z = new Tsil(  0,   0,   0,  0, MMZ, mu2);
  Wprotos[20] =  Wprot00000 = new Tsil(  0,   0,   0,   0,  0, mu2);
  Wprotos[21] =  WprotHW00  = new TsilSTU(MMH, MMW,   0,   0, mu2);
  Wprotos[22] =  WprotWH0H  = new TsilSTU(MMW, MMH,   0, MMH, mu2);
  Wprotos[23] =  WprotWZ0Z  = new TsilSTU(MMW, MMZ,   0, MMZ, mu2);

  // 
  Zprotos[0] = ZprotZHHZZ = new Tsil(MMZ, MMH, MMH, MMZ, MMZ, mu2);
  Zprotos[1] = ZprotZZHHH = new Tsil(MMZ, MMZ, MMH, MMH, MMH, mu2);
  Zprotos[2] = ZprotZWHWW = new Tsil(MMZ, MMW, MMH, MMW, MMW, mu2);
  Zprotos[3] = ZprottZtHt = new Tsil(MMt, MMZ, MMt, MMH, MMt, mu2);
  Zprotos[4] = ZprotWWWWH = new Tsil(MMW, MMW, MMW, MMW, MMH, mu2);
  Zprotos[5] = ZprotWWWWZ = new Tsil(MMW, MMW, MMW, MMW, MMZ, mu2);
  Zprotos[6] = ZprotWWWW0 = new Tsil(MMW, MMW, MMW, MMW,   0, mu2);
  Zprotos[7] = ZprotWtWt0 = new Tsil(MMW, MMt, MMW, MMt,   0, mu2);
  Zprotos[8] = ZprotW0W0t = new Tsil(MMW,   0, MMW,   0, MMt, mu2);
  Zprotos[9] = ZprotW0W00 = new Tsil(MMW,   0, MMW,   0,   0, mu2);
  Zprotos[10] = ZprotttttH = new Tsil(MMt, MMt, MMt, MMt, MMH, mu2);
  Zprotos[11] = ZprotttttZ = new Tsil(MMt, MMt, MMt, MMt, MMZ, mu2);
  Zprotos[12] = Zprottttt0 = new Tsil(MMt, MMt, MMt, MMt,   0, mu2);
  Zprotos[13] = Zprott0t0W = new Tsil(MMt,   0, MMt,   0, MMW, mu2);
  Zprotos[14] = Zprot0000Z = new Tsil(  0,   0,   0,   0, MMZ,
  mu2);
  Zprotos[15] = Zprot0000W = new Tsil(  0,   0,   0,   0, MMW, mu2);
  Zprotos[16] = Zprot00000 = new Tsil(  0,   0,   0,   0,   0,
  mu2);
  Zprotos[17] = ZprotWZWHW = new Tsil(MMW, MMZ, MMW, MMH, MMW, mu2);
  Zprotos[18] = ZprotHZ00  = new TsilSTU(MMH, MMZ,   0,  0, mu2);

  
  Timer t1;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(int i = 0 ; i < 24; i++)
    Wprotos[i]->evaluate(MMW);
  
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(int i = 0 ; i < 19; i++)
    Zprotos[i]->evaluate(MMZ);
  t1.elapsed();
  
}


