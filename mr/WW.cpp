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

#include <WW.hpp>
#include "timer.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mr
{
  WW<OS>::WW(OSinput sm, long double mu2_)
  {
    MMb = sm.MMb();
    MMW = sm.MMW();
    MMZ = sm.MMZ();
    MMH = sm.MMH();
    MMt = sm.MMt();
    mu2 = mu2_;

    init();
  }


  void WW<OS>::init()
  {

    CW = sqrt(MMW/MMZ);
    SW = sqrt(1-MMW/MMZ);

    protWHHWW = std::unique_ptr<Tsil>(new Tsil(MMW, MMH, MMH, MMW, MMW, mu2) );
    protWHZWW = std::unique_ptr<Tsil>(new Tsil(MMW, MMH, MMZ, MMW, MMW, mu2) );
    protWZZWW = std::unique_ptr<Tsil>(new Tsil(MMW, MMZ, MMZ, MMW, MMW, mu2) );
    protWWHHH = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMH, MMH, MMH, mu2) );
    protWWHZZ = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMH, MMZ, MMZ, mu2) );
    protWWZZH = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMZ, MMZ, MMH, mu2) );
    protWtZ00 = std::unique_ptr<Tsil>(new Tsil(MMW, MMt, MMZ,   0,   0, mu2) );
    protW0HWW = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMH, MMW, MMW, mu2) );
    protW0Htt = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMH, MMt, MMt, mu2) );
    protW0ZWW = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMZ, MMW, MMW, mu2) );
    protW0Ztt = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMZ, MMt, MMt, mu2) );
    protW0Z00 = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMZ,   0,   0, mu2) );
    prot0WW0W = std::unique_ptr<Tsil>(new Tsil(  0, MMW, MMW,   0, MMW, mu2) );
    prot0Wt0t = std::unique_ptr<Tsil>(new Tsil(  0, MMW, MMt,   0, MMt, mu2) );
    prot0W0Z0 = std::unique_ptr<Tsil>(new Tsil(  0, MMW,   0, MMZ,   0, mu2) );
    prot00Wt0 = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMW, MMt,   0, mu2) );
    prot00W00 = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMW,   0,   0, mu2) );
    prot00ttZ = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMt, MMt, MMZ, mu2) );
    prot00tt0 = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMt, MMt,   0, mu2) );
    prot0000Z = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,  0, MMZ, mu2) );
    prot00000 = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0,  0, mu2) );
    protWH0H = std::unique_ptr<TsilSTU>(new TsilSTU(MMW, MMH,   0, MMH, mu2) );
    protWZ0Z = std::unique_ptr<TsilSTU>(new TsilSTU(MMW, MMZ,   0, MMZ, mu2) );
    protHW00 = std::unique_ptr<TsilSTU>(new TsilSTU(MMH, MMW,   0,   0, mu2) );

    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
      protWHHWW->evaluate(MMW);
#pragma omp section
      protWHZWW->evaluate(MMW);
#pragma omp section
      protWZZWW->evaluate(MMW);
#pragma omp section
      protWWHHH->evaluate(MMW);
#pragma omp section
      protWWHZZ->evaluate(MMW);
#pragma omp section
      protWWZZH->evaluate(MMW);
#pragma omp section
      protWtZ00->evaluate(MMW);
#pragma omp section
      protW0HWW->evaluate(MMW);
#pragma omp section
      protW0Htt->evaluate(MMW);
#pragma omp section
      protW0ZWW->evaluate(MMW);
#pragma omp section
      protW0Ztt->evaluate(MMW);
#pragma omp section
      protW0Z00->evaluate(MMW);
#pragma omp section
      prot0WW0W->evaluate(MMW);
#pragma omp section
      prot0Wt0t->evaluate(MMW);
#pragma omp section
      prot0W0Z0->evaluate(MMW);
#pragma omp section
      prot00Wt0->evaluate(MMW);
#pragma omp section
      prot00W00->evaluate(MMW);
#pragma omp section
      prot00ttZ->evaluate(MMW);
#pragma omp section
      prot00tt0->evaluate(MMW);
#pragma omp section
      prot0000Z->evaluate(MMW);
#pragma omp section
      prot00000->evaluate(MMW);
#pragma omp section
      protWH0H->evaluate(MMW);
#pragma omp section
      protWZ0Z->evaluate(MMW);
#pragma omp section
      protHW00->evaluate(MMW);
    }
#else
    protWHHWW->evaluate(MMW);
    protWHZWW->evaluate(MMW);
    protWZZWW->evaluate(MMW);
    protWWHHH->evaluate(MMW);
    protWWHZZ->evaluate(MMW);
    protWWZZH->evaluate(MMW);
    protWtZ00->evaluate(MMW);
    protW0HWW->evaluate(MMW);
    protW0Htt->evaluate(MMW);
    protW0ZWW->evaluate(MMW);
    protW0Ztt->evaluate(MMW);
    protW0Z00->evaluate(MMW);
    prot0WW0W->evaluate(MMW);
    prot0Wt0t->evaluate(MMW);
    prot0W0Z0->evaluate(MMW);
    prot00Wt0->evaluate(MMW);
    prot00W00->evaluate(MMW);
    prot00ttZ->evaluate(MMW);
    prot00tt0->evaluate(MMW);
    prot0000Z->evaluate(MMW);
    prot00000->evaluate(MMW);
    protWH0H->evaluate(MMW);
    protWZ0Z->evaluate(MMW);
    protHW00->evaluate(MMW);
#endif
    t.elapsed();
  }


  WW<MS>::WW(MSinput sm, long double mu2_)
  {
    mmb = sm.mmb();
    mmW = sm.mmW();
    mmZ = sm.mmZ();
    mmH = sm.mmH();
    mmt = sm.mmt();
    mu2 = mu2_;

    init();
  }


  void WW<MS>::init()
  {

    c = sqrt(mmW/mmZ);
    s = sqrt(1-mmW/mmZ);

    protWHHWW = std::unique_ptr<Tsil>(new Tsil(mmW, mmH, mmH, mmW, mmW, mu2) );
    protWHZWW = std::unique_ptr<Tsil>(new Tsil(mmW, mmH, mmZ, mmW, mmW, mu2) );
    protWZZWW = std::unique_ptr<Tsil>(new Tsil(mmW, mmZ, mmZ, mmW, mmW, mu2) );
    protWWHHH = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmH, mmH, mmH, mu2) );
    protWWHZZ = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmH, mmZ, mmZ, mu2) );
    protWWZZH = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmZ, mmZ, mmH, mu2) );
    protWtZ00 = std::unique_ptr<Tsil>(new Tsil(mmW, mmt, mmZ,   0,   0, mu2) );
    protW0HWW = std::unique_ptr<Tsil>(new Tsil(mmW,   0, mmH, mmW, mmW, mu2) );
    protW0Htt = std::unique_ptr<Tsil>(new Tsil(mmW,   0, mmH, mmt, mmt, mu2) );
    protW0ZWW = std::unique_ptr<Tsil>(new Tsil(mmW,   0, mmZ, mmW, mmW, mu2) );
    protW0Ztt = std::unique_ptr<Tsil>(new Tsil(mmW,   0, mmZ, mmt, mmt, mu2) );
    protW0Z00 = std::unique_ptr<Tsil>(new Tsil(mmW,   0, mmZ,   0,   0, mu2) );
    prot0WW0W = std::unique_ptr<Tsil>(new Tsil(  0, mmW, mmW,   0, mmW, mu2) );
    prot0Wt0t = std::unique_ptr<Tsil>(new Tsil(  0, mmW, mmt,   0, mmt, mu2) );
    prot0W0Z0 = std::unique_ptr<Tsil>(new Tsil(  0, mmW,   0, mmZ,   0, mu2) );
    prot00Wt0 = std::unique_ptr<Tsil>(new Tsil(  0,   0, mmW, mmt,   0, mu2) );
    prot00W00 = std::unique_ptr<Tsil>(new Tsil(  0,   0, mmW,   0,   0, mu2) );
    prot00ttZ = std::unique_ptr<Tsil>(new Tsil(  0,   0, mmt, mmt, mmZ, mu2) );
    prot00tt0 = std::unique_ptr<Tsil>(new Tsil(  0,   0, mmt, mmt,   0, mu2) );
    prot0000Z = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,  0, mmZ, mu2) );
    prot00000 = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0,  0, mu2) );
    protWH0H = std::unique_ptr<TsilSTU>(new TsilSTU(mmW, mmH,   0, mmH, mu2) );
    protWZ0Z = std::unique_ptr<TsilSTU>(new TsilSTU(mmW, mmZ,   0, mmZ, mu2) );
    protHW00 = std::unique_ptr<TsilSTU>(new TsilSTU(mmH, mmW,   0,   0, mu2) );


    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
    protWHHWW->evaluate(mmW);
#pragma omp section
    protWHZWW->evaluate(mmW);
#pragma omp section
    protWZZWW->evaluate(mmW);
#pragma omp section
    protWWHHH->evaluate(mmW);
#pragma omp section
    protWWHZZ->evaluate(mmW);
#pragma omp section
    protWWZZH->evaluate(mmW);
#pragma omp section
    protWtZ00->evaluate(mmW);
#pragma omp section
    protW0HWW->evaluate(mmW);
#pragma omp section
    protW0Htt->evaluate(mmW);
#pragma omp section
    protW0ZWW->evaluate(mmW);
#pragma omp section
    protW0Ztt->evaluate(mmW);
#pragma omp section
    protW0Z00->evaluate(mmW);
#pragma omp section
    prot0WW0W->evaluate(mmW);
#pragma omp section
    prot0Wt0t->evaluate(mmW);
#pragma omp section
    prot0W0Z0->evaluate(mmW);
#pragma omp section
    prot00Wt0->evaluate(mmW);
#pragma omp section
    prot00W00->evaluate(mmW);
#pragma omp section
    prot00ttZ->evaluate(mmW);
#pragma omp section
    prot00tt0->evaluate(mmW);
#pragma omp section
    prot0000Z->evaluate(mmW);
#pragma omp section
    prot00000->evaluate(mmW);
#pragma omp section
    protWH0H->evaluate(mmW);
#pragma omp section
    protWZ0Z->evaluate(mmW);
#pragma omp section
    protHW00->evaluate(mmW);
    }
#else
    protWHHWW->evaluate(mmW);
    protWHZWW->evaluate(mmW);
    protWZZWW->evaluate(mmW);
    protWWHHH->evaluate(mmW);
    protWWHZZ->evaluate(mmW);
    protWWZZH->evaluate(mmW);
    protWtZ00->evaluate(mmW);
    protW0HWW->evaluate(mmW);
    protW0Htt->evaluate(mmW);
    protW0ZWW->evaluate(mmW);
    protW0Ztt->evaluate(mmW);
    protW0Z00->evaluate(mmW);
    prot0WW0W->evaluate(mmW);
    prot0Wt0t->evaluate(mmW);
    prot0W0Z0->evaluate(mmW);
    prot00Wt0->evaluate(mmW);
    prot00W00->evaluate(mmW);
    prot00ttZ->evaluate(mmW);
    prot00tt0->evaluate(mmW);
    prot0000Z->evaluate(mmW);
    prot00000->evaluate(mmW);
    protWH0H->evaluate(mmW);
    protWZ0Z->evaluate(mmW);
    protHW00->evaluate(mmW);
#endif
    t.elapsed();
  }
} // namespace mr

