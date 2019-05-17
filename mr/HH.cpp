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

#include <HH.hpp>
#include "timer.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mr
{
  HH<OS>::HH(OSinput sm, long double mu2_)
  {
    MMb = sm.MMb();
    MMW = sm.MMW();
    MMZ = sm.MMZ();
    MMH = sm.MMH();
    MMt = sm.MMt();
    mu2 = mu2_;

    init();
  }


  void HH<OS>::init()
  {
  
    CW = sqrt(MMW/MMZ);
    SW = sqrt(1-MMW/MMZ);
  
    protHHHHH = std::unique_ptr<Tsil>(new Tsil(MMH, MMH, MMH, MMH, MMH, mu2) );
    protHZHZZ = std::unique_ptr<Tsil>(new Tsil(MMH, MMZ, MMH, MMZ, MMZ, mu2) );
    protHWHWW = std::unique_ptr<Tsil>(new Tsil(MMH, MMW, MMH, MMW, MMW, mu2) );
    protHtHtt = std::unique_ptr<Tsil>(new Tsil(MMH, MMt, MMH, MMt, MMt, mu2) );
    protZZZZH = std::unique_ptr<Tsil>(new Tsil(MMZ, MMZ, MMZ, MMZ, MMH, mu2) );
    protZWZWW = std::unique_ptr<Tsil>(new Tsil(MMZ, MMW, MMZ, MMW, MMW, mu2) );
    protZtZtt = std::unique_ptr<Tsil>(new Tsil(MMZ, MMt, MMZ, MMt, MMt, mu2) );
    protWWWWH = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW, MMH, mu2) );
    protWWWWZ = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW, MMZ, mu2) );
    protWWWW0 = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW,   0, mu2) );
    protWtWt0 = std::unique_ptr<Tsil>(new Tsil(MMW, MMt, MMW, MMt,   0, mu2) );
    protttttH = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt, MMH, mu2) );
    protttttZ = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt, MMZ, mu2) );
    prottttt0 = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt,   0, mu2) );
    protZZ00 = std::unique_ptr<TsilSTU>(new TsilSTU(MMZ, MMZ,    0,   0, mu2) );
    protWW00 = std::unique_ptr<TsilSTU>(new TsilSTU(MMW, MMW,    0,   0, mu2) );

    prot0H0H0 = std::unique_ptr<Tsil>(new Tsil(  0, MMH,   0, MMH,   0, mu2) );
    prot0t0tt = std::unique_ptr<Tsil>(new Tsil(  0, MMt,   0, MMt, MMt, mu2) );
    prot0t0t0 = std::unique_ptr<Tsil>(new Tsil(  0, MMt,   0, MMt,   0, mu2) );
    prot0000H = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0, MMH, mu2) );

    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
    protHHHHH->evaluate(MMH);
#pragma omp section
    protHZHZZ->evaluate(MMH);
#pragma omp section
    protHWHWW->evaluate(MMH);
#pragma omp section
    protHtHtt->evaluate(MMH);
#pragma omp section
    protZZZZH->evaluate(MMH);
#pragma omp section
    protZWZWW->evaluate(MMH);
#pragma omp section
    protZtZtt->evaluate(MMH);
#pragma omp section
    protWWWWH->evaluate(MMH);
#pragma omp section
    protWWWWZ->evaluate(MMH);
#pragma omp section
    protWWWW0->evaluate(MMH);
#pragma omp section
    protWtWt0->evaluate(MMH);
#pragma omp section
    protttttH->evaluate(MMH);
#pragma omp section
    protttttZ->evaluate(MMH);
#pragma omp section
    prottttt0->evaluate(MMH);
#pragma omp section
    protZZ00 ->evaluate(MMH);
#pragma omp section
    protWW00 ->evaluate(MMH);

#pragma omp section
    prot0H0H0->evaluate(MMH);
#pragma omp section
    prot0t0tt->evaluate(MMH);
#pragma omp section
    prot0t0t0->evaluate(MMH);
#pragma omp section
    prot0000H->evaluate(MMH);
    }
#else
    protHHHHH->evaluate(MMH);
    protHZHZZ->evaluate(MMH);
    protHWHWW->evaluate(MMH);
    protHtHtt->evaluate(MMH);
    protZZZZH->evaluate(MMH);
    protZWZWW->evaluate(MMH);
    protZtZtt->evaluate(MMH);
    protWWWWH->evaluate(MMH);
    protWWWWZ->evaluate(MMH);
    protWWWW0->evaluate(MMH);
    protWtWt0->evaluate(MMH);
    protttttH->evaluate(MMH);
    protttttZ->evaluate(MMH);
    prottttt0->evaluate(MMH);
    protZZ00 ->evaluate(MMH);
    protWW00 ->evaluate(MMH);

    prot0H0H0->evaluate(MMH);
    prot0t0tt->evaluate(MMH);
    prot0t0t0->evaluate(MMH);
    prot0000H->evaluate(MMH);
#endif
    t.elapsed();
  }



  HH<MS>::HH(MSinput sm, long double mu2_)
  {
    mmb = sm.mmb();
    mmW = sm.mmW();
    mmZ = sm.mmZ();
    mmH = sm.mmH();
    mmt = sm.mmt();
    mu2 = mu2_;

    init();
  }


  void HH<MS>::init()
  {
  
    c = sqrt(mmW/mmZ);
    s = sqrt(1-mmW/mmZ);
  
    protHHHHH = std::unique_ptr<Tsil>(new Tsil(mmH, mmH, mmH, mmH, mmH, mu2) );
    protHZHZZ = std::unique_ptr<Tsil>(new Tsil(mmH, mmZ, mmH, mmZ, mmZ, mu2) );
    protHWHWW = std::unique_ptr<Tsil>(new Tsil(mmH, mmW, mmH, mmW, mmW, mu2) );
    protHtHtt = std::unique_ptr<Tsil>(new Tsil(mmH, mmt, mmH, mmt, mmt, mu2) );
    protZZZZH = std::unique_ptr<Tsil>(new Tsil(mmZ, mmZ, mmZ, mmZ, mmH, mu2) );
    protZWZWW = std::unique_ptr<Tsil>(new Tsil(mmZ, mmW, mmZ, mmW, mmW, mu2) );
    protZtZtt = std::unique_ptr<Tsil>(new Tsil(mmZ, mmt, mmZ, mmt, mmt, mu2) );
    protWWWWH = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmW, mmW, mmH, mu2) );
    protWWWWZ = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmW, mmW, mmZ, mu2) );
    protWWWW0 = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmW, mmW,   0, mu2) );
    protWtWt0 = std::unique_ptr<Tsil>(new Tsil(mmW, mmt, mmW, mmt,   0, mu2) );
    protttttH = std::unique_ptr<Tsil>(new Tsil(mmt, mmt, mmt, mmt, mmH, mu2) );
    protttttZ = std::unique_ptr<Tsil>(new Tsil(mmt, mmt, mmt, mmt, mmZ, mu2) );
    prottttt0 = std::unique_ptr<Tsil>(new Tsil(mmt, mmt, mmt, mmt,   0, mu2) );
    protZZ00 = std::unique_ptr<TsilSTU>(new TsilSTU(mmZ, mmZ,    0,   0, mu2) );
    protWW00 = std::unique_ptr<TsilSTU>(new TsilSTU(mmW, mmW,    0,   0, mu2) );

    prot0H0H0 = std::unique_ptr<Tsil>(new Tsil(  0, mmH,   0, mmH,   0, mu2) );
    prot0t0tt = std::unique_ptr<Tsil>(new Tsil(  0, mmt,   0, mmt, mmt, mu2) );
    prot0t0t0 = std::unique_ptr<Tsil>(new Tsil(  0, mmt,   0, mmt,   0, mu2) );
    prot0000H = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0, mmH, mu2) );



    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
    protHHHHH->evaluate(mmH);
#pragma omp section
    protHZHZZ->evaluate(mmH);
#pragma omp section
    protHWHWW->evaluate(mmH);
#pragma omp section
    protHtHtt->evaluate(mmH);
#pragma omp section
    protZZZZH->evaluate(mmH);
#pragma omp section
    protZWZWW->evaluate(mmH);
#pragma omp section
    protZtZtt->evaluate(mmH);
#pragma omp section
    protWWWWH->evaluate(mmH);
#pragma omp section
    protWWWWZ->evaluate(mmH);
#pragma omp section
    protWWWW0->evaluate(mmH);
#pragma omp section
    protWtWt0->evaluate(mmH);
#pragma omp section
    protttttH->evaluate(mmH);
#pragma omp section
    protttttZ->evaluate(mmH);
#pragma omp section
    prottttt0->evaluate(mmH);
#pragma omp section
    protZZ00->evaluate(mmH);
#pragma omp section
    protWW00->evaluate(mmH);

#pragma omp section
    prot0H0H0->evaluate(mmH);
#pragma omp section
    prot0t0tt->evaluate(mmH);
#pragma omp section
    prot0t0t0->evaluate(mmH);
#pragma omp section
    prot0000H->evaluate(mmH);
    }
#else
    protHHHHH->evaluate(mmH);
    protHZHZZ->evaluate(mmH);
    protHWHWW->evaluate(mmH);
    protHtHtt->evaluate(mmH);
    protZZZZH->evaluate(mmH);
    protZWZWW->evaluate(mmH);
    protZtZtt->evaluate(mmH);
    protWWWWH->evaluate(mmH);
    protWWWWZ->evaluate(mmH);
    protWWWW0->evaluate(mmH);
    protWtWt0->evaluate(mmH);
    protttttH->evaluate(mmH);
    protttttZ->evaluate(mmH);
    prottttt0->evaluate(mmH);
    protZZ00->evaluate(mmH);
    protWW00->evaluate(mmH);

    prot0H0H0->evaluate(mmH);
    prot0t0tt->evaluate(mmH);
    prot0t0t0->evaluate(mmH);
    prot0000H->evaluate(mmH);
#endif
    t.elapsed();
  }
} // namespace mr
