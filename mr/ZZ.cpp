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

#include <ZZ.hpp>
#include "timer.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mr
{
  ZZ<OS>::ZZ(OSinput sm, long double mu2_)
  {
    MMb = sm.MMb();
    MMW = sm.MMW();
    MMZ = sm.MMZ();
    MMH = sm.MMH();
    MMt = sm.MMt();
    mu2 = mu2_;

    init();
  }

  void ZZ<OS>::init()
  {

    CW = sqrt(MMW/MMZ);
    SW = sqrt(1-MMW/MMZ);

    protZHHZZ = std::unique_ptr<Tsil>(new Tsil(MMZ, MMH, MMH, MMZ, MMZ, mu2) );
    protZZHHH = std::unique_ptr<Tsil>(new Tsil(MMZ, MMZ, MMH, MMH, MMH, mu2) );
    protZWHWW = std::unique_ptr<Tsil>(new Tsil(MMZ, MMW, MMH, MMW, MMW, mu2) );
    prottZtHt = std::unique_ptr<Tsil>(new Tsil(MMt, MMZ, MMt, MMH, MMt, mu2) );
    protWWWWH = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW, MMH, mu2) );
    protWWWWZ = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW, MMZ, mu2) );
    protWWWW0 = std::unique_ptr<Tsil>(new Tsil(MMW, MMW, MMW, MMW,   0, mu2) );
    protWtWt0 = std::unique_ptr<Tsil>(new Tsil(MMW, MMt, MMW, MMt,   0, mu2) );
    protW0W0t = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMW,   0, MMt, mu2) );
    protW0W00 = std::unique_ptr<Tsil>(new Tsil(MMW,   0, MMW,   0,   0, mu2) );
    protttttH = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt, MMH, mu2) );
    protttttZ = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt, MMZ, mu2) );
    prottttt0 = std::unique_ptr<Tsil>(new Tsil(MMt, MMt, MMt, MMt,   0, mu2) );
    prott0t0W = std::unique_ptr<Tsil>(new Tsil(MMt,   0, MMt,   0, MMW, mu2) );
    prot0000Z = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0, MMZ, mu2) );
    prot0000W = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0, MMW, mu2) );
    prot00000 = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0,   0, mu2) );
    protWZWHW = std::unique_ptr<Tsil>(new Tsil(MMW, MMZ, MMW, MMH, MMW, mu2) );
    protHZ00  = std::unique_ptr<TsilSTU>(new TsilSTU(MMH, MMZ,   0,  0, mu2) );

    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
    protZHHZZ->evaluate(MMZ);
#pragma omp section
    protZZHHH->evaluate(MMZ);
#pragma omp section
    protZWHWW->evaluate(MMZ);
#pragma omp section
    prottZtHt->evaluate(MMZ);
#pragma omp section
    protWWWWH->evaluate(MMZ);
#pragma omp section
    protWWWWZ->evaluate(MMZ);
#pragma omp section
    protWWWW0->evaluate(MMZ);
#pragma omp section
    protWtWt0->evaluate(MMZ);
#pragma omp section
    protW0W0t->evaluate(MMZ);
#pragma omp section
    protW0W00->evaluate(MMZ);
#pragma omp section
    protttttH->evaluate(MMZ);
#pragma omp section
    protttttZ->evaluate(MMZ);
#pragma omp section
    prottttt0->evaluate(MMZ);
#pragma omp section
    prott0t0W->evaluate(MMZ);
#pragma omp section
    prot0000Z->evaluate(MMZ);
#pragma omp section
    prot0000W->evaluate(MMZ);
#pragma omp section
    prot00000->evaluate(MMZ);
#pragma omp section
    protWZWHW->evaluate(MMZ);
#pragma omp section
    protHZ00 ->evaluate(MMZ);
    }
#else
    protZHHZZ->evaluate(MMZ);
    protZZHHH->evaluate(MMZ);
    protZWHWW->evaluate(MMZ);
    prottZtHt->evaluate(MMZ);
    protWWWWH->evaluate(MMZ);
    protWWWWZ->evaluate(MMZ);
    protWWWW0->evaluate(MMZ);
    protWtWt0->evaluate(MMZ);
    protW0W0t->evaluate(MMZ);
    protW0W00->evaluate(MMZ);
    protttttH->evaluate(MMZ);
    protttttZ->evaluate(MMZ);
    prottttt0->evaluate(MMZ);
    prott0t0W->evaluate(MMZ);
    prot0000Z->evaluate(MMZ);
    prot0000W->evaluate(MMZ);
    prot00000->evaluate(MMZ);
    protWZWHW->evaluate(MMZ);
    protHZ00 ->evaluate(MMZ);
#endif
    t.elapsed();
  }


  ZZ<MS>::ZZ(MSinput sm, long double mu2_)
  {
    mmb = sm.mmb();
    mmW = sm.mmW();
    mmZ = sm.mmZ();
    mmH = sm.mmH();
    mmt = sm.mmt();
    mu2 = mu2_;

    init();
  }


  void ZZ<MS>::init()
  {

    c = sqrt(mmW/mmZ);
    s = sqrt(1-mmW/mmZ);

    protZHHZZ = std::unique_ptr<Tsil>(new Tsil(mmZ, mmH, mmH, mmZ, mmZ, mu2) );
    protZZHHH = std::unique_ptr<Tsil>(new Tsil(mmZ, mmZ, mmH, mmH, mmH, mu2) );
    protZWHWW = std::unique_ptr<Tsil>(new Tsil(mmZ, mmW, mmH, mmW, mmW, mu2) );
    prottZtHt = std::unique_ptr<Tsil>(new Tsil(mmt, mmZ, mmt, mmH, mmt, mu2) );
    protWWWWH = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmW, mmW, mmH, mu2) );
    protWWWWZ = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmW, mmW, mmZ, mu2) );
    protWWWW0 = std::unique_ptr<Tsil>(new Tsil(mmW, mmW, mmW, mmW,   0, mu2) );
    protWtWt0 = std::unique_ptr<Tsil>(new Tsil(mmW, mmt, mmW, mmt,   0, mu2) );
    protW0W0t = std::unique_ptr<Tsil>(new Tsil(mmW,   0, mmW,   0, mmt, mu2) );
    protW0W00 = std::unique_ptr<Tsil>(new Tsil(mmW,   0, mmW,   0,   0, mu2) );
    protttttH = std::unique_ptr<Tsil>(new Tsil(mmt, mmt, mmt, mmt, mmH, mu2) );
    protttttZ = std::unique_ptr<Tsil>(new Tsil(mmt, mmt, mmt, mmt, mmZ, mu2) );
    prottttt0 = std::unique_ptr<Tsil>(new Tsil(mmt, mmt, mmt, mmt,   0, mu2) );
    prott0t0W = std::unique_ptr<Tsil>(new Tsil(mmt,   0, mmt,   0, mmW, mu2) );
    prot0000Z = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0, mmZ, mu2) );
    prot0000W = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0, mmW, mu2) );
    prot00000 = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0,   0,   0, mu2) );
    protWZWHW = std::unique_ptr<Tsil>(new Tsil(mmW, mmZ, mmW, mmH, mmW, mu2) );
    protHZ00  = std::unique_ptr<TsilSTU>(new TsilSTU(mmH, mmZ,   0,  0, mu2) );

    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
    protZHHZZ->evaluate(mmZ);
#pragma omp section
    protZZHHH->evaluate(mmZ);
#pragma omp section
    protZWHWW->evaluate(mmZ);
#pragma omp section
    prottZtHt->evaluate(mmZ);
#pragma omp section
    protWWWWH->evaluate(mmZ);
#pragma omp section
    protWWWWZ->evaluate(mmZ);
#pragma omp section
    protWWWW0->evaluate(mmZ);
#pragma omp section
    protWtWt0->evaluate(mmZ);
#pragma omp section
    protW0W0t->evaluate(mmZ);
#pragma omp section
    protW0W00->evaluate(mmZ);
#pragma omp section
    protttttH->evaluate(mmZ);
#pragma omp section
    protttttZ->evaluate(mmZ);
#pragma omp section
    prottttt0->evaluate(mmZ);
#pragma omp section
    prott0t0W->evaluate(mmZ);
#pragma omp section
    prot0000Z->evaluate(mmZ);
#pragma omp section
    prot0000W->evaluate(mmZ);
#pragma omp section
    prot00000->evaluate(mmZ);
#pragma omp section
    protWZWHW->evaluate(mmZ);
#pragma omp section
    protHZ00 ->evaluate(mmZ);
    }
#else
    protZHHZZ->evaluate(mmZ);
    protZZHHH->evaluate(mmZ);
    protZWHWW->evaluate(mmZ);
    prottZtHt->evaluate(mmZ);
    protWWWWH->evaluate(mmZ);
    protWWWWZ->evaluate(mmZ);
    protWWWW0->evaluate(mmZ);
    protWtWt0->evaluate(mmZ);
    protW0W0t->evaluate(mmZ);
    protW0W00->evaluate(mmZ);
    protttttH->evaluate(mmZ);
    protttttZ->evaluate(mmZ);
    prottttt0->evaluate(mmZ);
    prott0t0W->evaluate(mmZ);
    prot0000Z->evaluate(mmZ);
    prot0000W->evaluate(mmZ);
    prot00000->evaluate(mmZ);
    protWZWHW->evaluate(mmZ);
    protHZ00 ->evaluate(mmZ);
#endif
    t.elapsed();
  }
} // namespace mr
