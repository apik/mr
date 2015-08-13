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

#include <tt.hpp>
#include "timer.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mr
{
  tt<OS>::tt(OSinput sm, long double mu2_)
  {
    MMb = sm.MMb();
    MMW = sm.MMW();
    MMZ = sm.MMZ();
    MMH = sm.MMH();
    MMt = sm.MMt();
    mu2 = mu2_;

    init();
  }


  void tt<OS>::init()
  {
  
    CW = sqrt(MMW/MMZ);
    SW = sqrt(1-MMW/MMZ);
  
    protos[0]  = protWt000 = new Tsil(MMW, MMt,   0,   0,   0, mu2);
    protos[1]  = prot0ttHt = new Tsil(0  , MMt, MMt, MMH, MMt, mu2);
    protos[2]  = prot0ttZt = new Tsil(0  , MMt, MMt, MMZ, MMt, mu2);
    protos[3]  = prot0tt0t = new Tsil(0  , MMt, MMt,   0, MMt, mu2);
    protos[4]  = prottH0H =  new TsilSTU ( MMt, MMH,   0, MMH, mu2);
    protos[5]  = prottZ0Z =  new TsilSTU ( MMt, MMZ,   0, MMZ, mu2);
  
    protos[6]  = protHHttH = new Tsil(MMH, MMH, MMt, MMt, MMH, mu2);
    protos[7]  = protHZttZ = new Tsil(MMH, MMZ, MMt, MMt, MMZ, mu2);
    protos[8]  = protHWt0W = new Tsil(MMH, MMW, MMt,   0, MMW, mu2);
    protos[9]  = protHttHt = new Tsil(MMH, MMt, MMt, MMH, MMt, mu2);
    protos[10] = protHttZt = new Tsil(MMH, MMt, MMt, MMZ, MMt, mu2);
    protos[11] = protZZttH = new Tsil(MMZ, MMZ, MMt, MMt, MMH, mu2);
    protos[12] = protZWt0W = new Tsil(MMZ, MMW, MMt,   0, MMW, mu2);
    protos[13] = protZttZt = new Tsil(MMZ, MMt, MMt, MMZ, MMt, mu2);
    protos[14] = protZ0tW0 = new Tsil(MMZ,   0, MMt, MMW,   0, mu2);
    protos[15] = protWW00Z = new Tsil(MMW, MMW,   0,   0, MMZ, mu2);
    protos[16] = protW00tW = new Tsil(MMW,   0,   0, MMt, MMW, mu2);
    protos[17] = prot00WW0 = new Tsil(0  ,   0, MMW, MMW,   0, mu2);
    protos[18] = prot000   = new TsilST(0,             0,   0, mu2);
    protos[19] = prot0W00  = new TsilSTU(0,     MMW,   0,   0, mu2);

    protos[20] = protH0tt0 = new Tsil(MMH,   0, MMt, MMt,   0, mu2);
    protos[21] = protH0t00 = new Tsil(MMH,   0, MMt,   0,   0, mu2);
    protos[22] = prot0Htt0 = new Tsil(  0, MMH, MMt, MMt,   0, mu2);
    protos[23] = prot0H0t0 = new Tsil(  0, MMH,   0, MMt,   0, mu2);
    protos[24] = prot00ttH = new Tsil(  0,   0, MMt, MMt, MMH, mu2);
    protos[25] = protHtt0t = new Tsil(MMH, MMt, MMt,   0, MMt, mu2);
    protos[26] = prot00t00 = new Tsil(  0,   0, MMt,   0,   0, mu2);
    protos[27] = prot000t0 = new Tsil(  0,   0,   0, MMt,   0, mu2);

  
    Timer t;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(int i = 0 ; i < 28; i++)
      protos[i]->evaluate(MMt);
    t.elapsed();
  
  }


  tt<MS>::tt(MSinput sm, long double mu2_)
  {
    mmb = sm.mmb();
    mmW = sm.mmW();
    mmZ = sm.mmZ();
    mmH = sm.mmH();
    mmt = sm.mmt();
    mu2 = mu2_;

    init();
  }


  void tt<MS>::init()
  {
  
    c = sqrt(mmW/mmZ);
    s = sqrt(1-mmW/mmZ);
  
    protos[0]  = protWt000 = new Tsil(mmW, mmt,   0,   0,   0, mu2);
    protos[1]  = prot0ttHt = new Tsil(0  , mmt, mmt, mmH, mmt, mu2);
    protos[2]  = prot0ttZt = new Tsil(0  , mmt, mmt, mmZ, mmt, mu2);
    protos[3]  = prot0tt0t = new Tsil(0  , mmt, mmt,   0, mmt, mu2);
    protos[4]  = prottH0H =  new TsilSTU ( mmt, mmH,   0, mmH, mu2);
    protos[5]  = prottZ0Z =  new TsilSTU ( mmt, mmZ,   0, mmZ, mu2);
  
    protos[6]  = protHHttH = new Tsil(mmH, mmH, mmt, mmt, mmH, mu2);
    protos[7]  = protHZttZ = new Tsil(mmH, mmZ, mmt, mmt, mmZ, mu2);
    protos[8]  = protHWt0W = new Tsil(mmH, mmW, mmt,   0, mmW, mu2);
    protos[9]  = protHttHt = new Tsil(mmH, mmt, mmt, mmH, mmt, mu2);
    protos[10] = protHttZt = new Tsil(mmH, mmt, mmt, mmZ, mmt, mu2);
    protos[11] = protZZttH = new Tsil(mmZ, mmZ, mmt, mmt, mmH, mu2);
    protos[12] = protZWt0W = new Tsil(mmZ, mmW, mmt,   0, mmW, mu2);
    protos[13] = protZttZt = new Tsil(mmZ, mmt, mmt, mmZ, mmt, mu2);
    protos[14] = protZ0tW0 = new Tsil(mmZ,   0, mmt, mmW,   0, mu2);
    protos[15] = protWW00Z = new Tsil(mmW, mmW,   0,   0, mmZ, mu2);
    protos[16] = protW00tW = new Tsil(mmW,   0,   0, mmt, mmW, mu2);
    protos[17] = prot00WW0 = new Tsil(0  ,   0, mmW, mmW,   0, mu2);
    protos[18] = prot000   = new TsilST(0,             0,   0, mu2);
    protos[19] = prot0W00  = new TsilSTU(0,     mmW,   0,   0, mu2);

    protos[20]  = protH0tt0 = new Tsil(mmH,   0, mmt, mmt,   0, mu2);
    protos[21]  = protH0t00 = new Tsil(mmH,   0, mmt,   0,   0, mu2);
    protos[22]  = prot0Htt0 = new Tsil(  0, mmH, mmt, mmt,   0, mu2);
    protos[23]  = prot0H0t0 = new Tsil(  0, mmH,   0, mmt,   0, mu2);
    protos[24]  = prot00ttH = new Tsil(  0,   0, mmt, mmt, mmH, mu2);
    protos[25]  = protHtt0t = new Tsil(mmH, mmt, mmt,   0, mmt, mu2);
    protos[26]  = prot00t00 = new Tsil(  0,   0, mmt,   0,   0, mu2);
    protos[27]  = prot000t0 = new Tsil(  0,   0,   0, mmt,   0, mu2);

  
    Timer t;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
    for(int i = 0 ; i < 28; i++)
      protos[i]->evaluate(mmt);
    t.elapsed();
  
  }
} // namespace mr
