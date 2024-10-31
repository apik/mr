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


namespace mr
{
  alphaGF::alphaGF(OSinput sm, double mu2_)
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

    WprotWHHWW = std::unique_ptr<Tsil>(new Tsil(MMW, MMH, MMH, MMW, MMW, mu2) );
    WprotWHZWW = std::unique_ptr<Tsil>(new Tsil(MMW, MMH, MMZ, MMW, MMW, mu2) );
    WprotWZZWW = std::unique_ptr<Tsil>(new Tsil(MMW, MMZ, MMZ, MMW, MMW, mu2) );
    WprotWWHHH = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMH, MMH, MMH, mu2) );
    WprotWWHZZ = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMH, MMZ, MMZ, mu2) );
    WprotWWZZH = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMZ, MMZ, MMH, mu2) );
    WprotWtZ00 = std::unique_ptr<Tsil>(new Tsil(MMW, MMt, MMZ,   0,   0, mu2) );
    WprotW0HWW = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMH, MMW, MMW, mu2) );
    WprotW0Htt = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMH, MMt, MMt, mu2) );
    WprotW0ZWW = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMZ, MMW, MMW, mu2) );
    WprotW0Ztt = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMZ, MMt, MMt, mu2) );
    WprotW0Z00 = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMZ,   0,   0, mu2) );
    Wprot0WW0W = std::unique_ptr<Tsil>(new Tsil(  0, MMW, MMW,   0, MMW, mu2) );
    Wprot0Wt0t = std::unique_ptr<Tsil>(new Tsil(  0, MMW, MMt,   0, MMt, mu2) );
    Wprot0W0Z0 = std::unique_ptr<Tsil>(new Tsil(  0, MMW,   0, MMZ,   0, mu2) );
    Wprot00Wt0 = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMW, MMt,   0, mu2) );
    Wprot00W00 = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMW,   0,   0, mu2) );
    Wprot00ttZ = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMt, MMt, MMZ, mu2) );
    Wprot00tt0 = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMt, MMt,   0, mu2) );
    Wprot0000Z = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,  0, MMZ, mu2) );
    Wprot00000 = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0,  0, mu2) );
    WprotHW00  = std::unique_ptr<TsilSTU>(new TsilSTU(MMH, MMW,   0,   0, mu2) );
    WprotWH0H  = std::unique_ptr<TsilSTU>(new TsilSTU(MMW, MMH,   0, MMH, mu2) );
    WprotWZ0Z  = std::unique_ptr<TsilSTU>(new TsilSTU(MMW, MMZ,   0, MMZ, mu2) );

    // 
    ZprotZHHZZ = std::unique_ptr<Tsil>(new Tsil(MMZ, MMH, MMH, MMZ, MMZ, mu2) );
    ZprotZZHHH = std::unique_ptr<Tsil>(new Tsil(MMZ, MMZ, MMH, MMH, MMH, mu2) );
    ZprotZWHWW = std::unique_ptr<Tsil>(new Tsil(MMZ, MMW, MMH, MMW, MMW, mu2) );
    ZprottZtHt = std::unique_ptr<Tsil>(new Tsil(MMt, MMZ, MMt, MMH, MMt, mu2) );
    ZprotWWWWH = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW, MMH, mu2) );
    ZprotWWWWZ = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW, MMZ, mu2) );
    ZprotWWWW0 = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW,   0, mu2) );
    ZprotWtWt0 = std::unique_ptr<Tsil>(new Tsil(MMW, MMt, MMW, MMt,   0, mu2) );
    ZprotW0W0t = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMW,   0, MMt, mu2) );
    ZprotW0W00 = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMW,   0,   0, mu2) );
    ZprotttttH = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt, MMH, mu2) );
    ZprotttttZ = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt, MMZ, mu2) );
    Zprottttt0 = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt,   0, mu2) );
    Zprott0t0W = std::unique_ptr<Tsil>(new Tsil(MMt,   0, MMt,   0, MMW, mu2) );
    Zprot0000Z = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0, MMZ, mu2) );
    Zprot0000W = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0, MMW, mu2) );
    Zprot00000 = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0,   0, mu2) );
    ZprotWZWHW = std::unique_ptr<Tsil>(new Tsil(MMW, MMZ, MMW, MMH, MMW, mu2) );
    ZprotHZ00  = std::unique_ptr<TsilSTU>(new TsilSTU(MMH, MMZ,   0,  0, mu2) );


    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
    WprotWHHWW->evaluate(MMW);
#pragma omp section
    WprotWHZWW->evaluate(MMW);
#pragma omp section
    WprotWZZWW->evaluate(MMW);
#pragma omp section
    WprotWWHHH->evaluate(MMW);
#pragma omp section
    WprotWWHZZ->evaluate(MMW);
#pragma omp section
    WprotWWZZH->evaluate(MMW);
#pragma omp section
    WprotWtZ00->evaluate(MMW);
#pragma omp section
    WprotW0HWW->evaluate(MMW);
#pragma omp section
    WprotW0Htt->evaluate(MMW);
#pragma omp section
    WprotW0ZWW->evaluate(MMW);
#pragma omp section
    WprotW0Ztt->evaluate(MMW);
#pragma omp section
    WprotW0Z00->evaluate(MMW);
#pragma omp section
    Wprot0WW0W->evaluate(MMW);
#pragma omp section
    Wprot0Wt0t->evaluate(MMW);
#pragma omp section
    Wprot0W0Z0->evaluate(MMW);
#pragma omp section
    Wprot00Wt0->evaluate(MMW);
#pragma omp section
    Wprot00W00->evaluate(MMW);
#pragma omp section
    Wprot00ttZ->evaluate(MMW);
#pragma omp section
    Wprot00tt0->evaluate(MMW);
#pragma omp section
    Wprot0000Z->evaluate(MMW);
#pragma omp section
    Wprot00000->evaluate(MMW);
#pragma omp section
    WprotHW00 ->evaluate(MMW);
#pragma omp section
    WprotWH0H ->evaluate(MMW);
#pragma omp section
    WprotWZ0Z ->evaluate(MMW);

    // 
#pragma omp section
    ZprotZHHZZ->evaluate(MMZ);
#pragma omp section
    ZprotZZHHH->evaluate(MMZ);
#pragma omp section
    ZprotZWHWW->evaluate(MMZ);
#pragma omp section
    ZprottZtHt->evaluate(MMZ);
#pragma omp section
    ZprotWWWWH->evaluate(MMZ);
#pragma omp section
    ZprotWWWWZ->evaluate(MMZ);
#pragma omp section
    ZprotWWWW0->evaluate(MMZ);
#pragma omp section
    ZprotWtWt0->evaluate(MMZ);
#pragma omp section
    ZprotW0W0t->evaluate(MMZ);
#pragma omp section
    ZprotW0W00->evaluate(MMZ);
#pragma omp section
    ZprotttttH->evaluate(MMZ);
#pragma omp section
    ZprotttttZ->evaluate(MMZ);
#pragma omp section
    Zprottttt0->evaluate(MMZ);
#pragma omp section
    Zprott0t0W->evaluate(MMZ);
#pragma omp section
    Zprot0000Z->evaluate(MMZ);
#pragma omp section
    Zprot0000W->evaluate(MMZ);
#pragma omp section
    Zprot00000->evaluate(MMZ);
#pragma omp section
    ZprotWZWHW->evaluate(MMZ);
#pragma omp section
    ZprotHZ00 ->evaluate(MMZ);
    }
    #else
    WprotWHHWW->evaluate(MMW);
    WprotWHZWW->evaluate(MMW);
    WprotWZZWW->evaluate(MMW);
    WprotWWHHH->evaluate(MMW);
    WprotWWHZZ->evaluate(MMW);
    WprotWWZZH->evaluate(MMW);
    WprotWtZ00->evaluate(MMW);
    WprotW0HWW->evaluate(MMW);
    WprotW0Htt->evaluate(MMW);
    WprotW0ZWW->evaluate(MMW);
    WprotW0Ztt->evaluate(MMW);
    WprotW0Z00->evaluate(MMW);
    Wprot0WW0W->evaluate(MMW);
    Wprot0Wt0t->evaluate(MMW);
    Wprot0W0Z0->evaluate(MMW);
    Wprot00Wt0->evaluate(MMW);
    Wprot00W00->evaluate(MMW);
    Wprot00ttZ->evaluate(MMW);
    Wprot00tt0->evaluate(MMW);
    Wprot0000Z->evaluate(MMW);
    Wprot00000->evaluate(MMW);
    WprotHW00 ->evaluate(MMW);
    WprotWH0H ->evaluate(MMW);
    WprotWZ0Z ->evaluate(MMW);

    //
    ZprotZHHZZ->evaluate(MMZ);
    ZprotZZHHH->evaluate(MMZ);
    ZprotZWHWW->evaluate(MMZ);
    ZprottZtHt->evaluate(MMZ);
    ZprotWWWWH->evaluate(MMZ);
    ZprotWWWWZ->evaluate(MMZ);
    ZprotWWWW0->evaluate(MMZ);
    ZprotWtWt0->evaluate(MMZ);
    ZprotW0W0t->evaluate(MMZ);
    ZprotW0W00->evaluate(MMZ);
    ZprotttttH->evaluate(MMZ);
    ZprotttttZ->evaluate(MMZ);
    Zprottttt0->evaluate(MMZ);
    Zprott0t0W->evaluate(MMZ);
    Zprot0000Z->evaluate(MMZ);
    Zprot0000W->evaluate(MMZ);
    Zprot00000->evaluate(MMZ);
    ZprotWZWHW->evaluate(MMZ);
    ZprotHZ00 ->evaluate(MMZ);
#endif
    t.elapsed();
  }
} // namespace mr
