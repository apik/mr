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

// #include <omp.h>
#include <alphaGF.hpp>
#include "timer.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

// WW::WW(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
//   MMb(MMb_), MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
// {
//   init(MMB, MMW, MMZ, MMH, MMt, mu2);
// }

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

  protos[0] =  protWHHWW = new Tsil(MMW, MMH, MMH, MMW, MMW, mu2);
  protos[1] =  protWHZWW = new Tsil(MMW, MMH, MMZ, MMW, MMW, mu2);
  protos[2] =  protWZZWW = new Tsil(MMW, MMZ, MMZ, MMW, MMW, mu2);
  protos[3] =  protWWHHH = new Tsil(MMW, MMW, MMH, MMH, MMH, mu2);
  protos[4] =  protWWHZZ = new Tsil(MMW, MMW, MMH, MMZ, MMZ, mu2);
  protos[5] =  protWWZZH = new Tsil(MMW, MMW, MMZ, MMZ, MMH, mu2);
  protos[6] =  protWtZ00 = new Tsil(MMW, MMt, MMZ,   0,   0, mu2);
  protos[7] =  protW0HWW = new Tsil(MMW,   0, MMH, MMW, MMW, mu2);
  protos[8] =  protW0Htt = new Tsil(MMW,   0, MMH, MMt, MMt, mu2);
  protos[9] =  protW0ZWW = new Tsil(MMW,   0, MMZ, MMW, MMW, mu2);
  protos[10] =  protW0Ztt = new Tsil(MMW,   0, MMZ, MMt, MMt, mu2);
  protos[11] =  protW0Z00 = new Tsil(MMW,   0, MMZ,   0,   0, mu2);
  protos[12] =  prot0WW0W = new Tsil(  0, MMW, MMW,   0, MMW, mu2);
  protos[13] =  prot0Wt0t = new Tsil(  0, MMW, MMt,   0, MMt, mu2);
  protos[14] =  prot0W0Z0 = new Tsil(  0, MMW,   0, MMZ,   0, mu2);
  protos[15] =  prot00Wt0 = new Tsil(  0,   0, MMW, MMt,   0, mu2);
  protos[16] =  prot00W00 = new Tsil(  0,   0, MMW,   0,   0, mu2);
  protos[17] =  prot00ttZ = new Tsil(  0,   0, MMt, MMt, MMZ, mu2);
  protos[18] =  prot00tt0 = new Tsil(  0,   0, MMt, MMt,   0, mu2);
  protos[19] =  prot0000Z = new Tsil(  0,   0,   0,  0, MMZ, mu2);
  protos[20] =  prot00000 = new Tsil(  0,   0,   0,   0,  0, mu2);
  protos[21] =  protWH0H = new TsilSTU(MMW, MMH,   0, MMH, mu2);
  protos[22] =  protWZ0Z = new TsilSTU(MMW, MMZ,   0, MMZ, mu2);
  protos[23] =  protHW00 = new TsilSTU(MMH, MMW,   0,   0, mu2);
  // 
  protos[24] = protZHHZZ = new Tsil(MMZ, MMH, MMH, MMZ, MMZ, mu2);
  protos[25] = protZZHHH = new Tsil(MMZ, MMZ, MMH, MMH, MMH, mu2);
  protos[26] = protZWHWW = new Tsil(MMZ, MMW, MMH, MMW, MMW, mu2);
  protos[27] = prottZtHt = new Tsil(MMt, MMZ, MMt, MMH, MMt, mu2);
  protos[28] = protWWWWH = new Tsil(MMW, MMW, MMW, MMW, MMH, mu2);
  protos[29] = protWWWWZ = new Tsil(MMW, MMW, MMW, MMW, MMZ, mu2);
  protos[30] = protWWWW0 = new Tsil(MMW, MMW, MMW, MMW,   0, mu2);
  protos[31] = protWtWt0 = new Tsil(MMW, MMt, MMW, MMt,   0, mu2);
  protos[32] = protW0W0t = new Tsil(MMW,   0, MMW,   0, MMt, mu2);
  protos[33] = protW0W00 = new Tsil(MMW,   0, MMW,   0,   0, mu2);
  protos[34] = protttttH = new Tsil(MMt, MMt, MMt, MMt, MMH, mu2);
  protos[35] = protttttZ = new Tsil(MMt, MMt, MMt, MMt, MMZ, mu2);
  protos[36] = prottttt0 = new Tsil(MMt, MMt, MMt, MMt,   0, mu2);
  protos[37] = prott0t0W = new Tsil(MMt,   0, MMt,   0, MMW, mu2);

  protos[38] = prot0000W = new Tsil(  0,   0,   0,   0, MMW, mu2);

  protos[39] = protWZWHW = new Tsil(MMW, MMZ, MMW, MMH, MMW, mu2);
  protos[40] = protHZ00  = new TsilSTU(MMH, MMZ,   0,  0, mu2);

  
  Timer t1;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
  for(int i = 0 ; i < 41; i++)
    protos[i]->evaluate(MMW);
  t1.elapsed();
  
}


