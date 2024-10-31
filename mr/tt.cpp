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
  tt<OS>::tt(OSinput sm, double mu2_)
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
  
    protWt000 = std::unique_ptr<Tsil>(new Tsil(MMW, MMt,   0,   0,   0, mu2) );
    prot0ttHt = std::unique_ptr<Tsil>(new Tsil(0  , MMt, MMt, MMH, MMt, mu2) );
    prot0ttZt = std::unique_ptr<Tsil>(new Tsil(0  , MMt, MMt, MMZ, MMt, mu2) );
    prot0tt0t = std::unique_ptr<Tsil>(new Tsil(0  , MMt, MMt,   0, MMt, mu2) );
    prottH0H =  std::unique_ptr<TsilSTU>(new TsilSTU ( MMt, MMH,   0, MMH, mu2) );
    prottZ0Z =  std::unique_ptr<TsilSTU>(new TsilSTU ( MMt, MMZ,   0, MMZ, mu2) );
    
    protHHttH = std::unique_ptr<Tsil>(new Tsil(MMH, MMH, MMt, MMt, MMH, mu2) );
    protHZttZ = std::unique_ptr<Tsil>(new Tsil(MMH, MMZ, MMt, MMt, MMZ, mu2) );
    protHWt0W = std::unique_ptr<Tsil>(new Tsil(MMH, MMW, MMt,   0, MMW, mu2) );
    protHttHt = std::unique_ptr<Tsil>(new Tsil(MMH, MMt, MMt, MMH, MMt, mu2) );
    protHttZt = std::unique_ptr<Tsil>(new Tsil(MMH, MMt, MMt, MMZ, MMt, mu2) );
    protZZttH = std::unique_ptr<Tsil>(new Tsil(MMZ, MMZ, MMt, MMt, MMH, mu2) );
    protZWt0W = std::unique_ptr<Tsil>(new Tsil(MMZ, MMW, MMt,   0, MMW, mu2) );
    protZttZt = std::unique_ptr<Tsil>(new Tsil(MMZ, MMt, MMt, MMZ, MMt, mu2) );
    protZ0tW0 = std::unique_ptr<Tsil>(new Tsil(MMZ,   0, MMt, MMW,   0, mu2) );
    protWW00Z = std::unique_ptr<Tsil>(new Tsil(MMW, MMW,   0,   0, MMZ, mu2) );
    protW00tW = std::unique_ptr<Tsil>(new Tsil(MMW,   0,   0, MMt, MMW, mu2) );
    prot00WW0 = std::unique_ptr<Tsil>(new Tsil(0  ,   0, MMW, MMW,   0, mu2) );
    prot000   = std::unique_ptr<TsilST>(new TsilST(0,             0,   0, mu2) );
    prot0W00  = std::unique_ptr<TsilSTU>(new TsilSTU(0,     MMW,   0,   0, mu2) );
    
    protH0tt0 = std::unique_ptr<Tsil>(new Tsil(MMH,   0, MMt, MMt,   0, mu2) );
    protH0t00 = std::unique_ptr<Tsil>(new Tsil(MMH,   0, MMt,   0,   0, mu2) );
    prot0Htt0 = std::unique_ptr<Tsil>(new Tsil(  0, MMH, MMt, MMt,   0, mu2) );
    prot0H0t0 = std::unique_ptr<Tsil>(new Tsil(  0, MMH,   0, MMt,   0, mu2) );
    prot00ttH = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMt, MMt, MMH, mu2) );
    protHtt0t = std::unique_ptr<Tsil>(new Tsil(MMH, MMt, MMt,   0, MMt, mu2) );
    prot00t00 = std::unique_ptr<Tsil>(new Tsil(  0,   0, MMt,   0,   0, mu2) );
    prot000t0 = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0, MMt,   0, mu2) );


    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
    protWt000->evaluate(MMt);
#pragma omp section
    prot0ttHt->evaluate(MMt);
#pragma omp section
    prot0ttZt->evaluate(MMt);
#pragma omp section
    prot0tt0t->evaluate(MMt);
#pragma omp section
    prottH0H->evaluate(MMt);
#pragma omp section
    prottZ0Z->evaluate(MMt);

#pragma omp section
    protHHttH->evaluate(MMt);
#pragma omp section
    protHZttZ->evaluate(MMt);
#pragma omp section
    protHWt0W->evaluate(MMt);
#pragma omp section
    protHttHt->evaluate(MMt);
#pragma omp section
    protHttZt->evaluate(MMt);
#pragma omp section
    protZZttH->evaluate(MMt);
#pragma omp section
    protZWt0W->evaluate(MMt);
#pragma omp section
    protZttZt->evaluate(MMt);
#pragma omp section
    protZ0tW0->evaluate(MMt);
#pragma omp section
    protWW00Z->evaluate(MMt);
#pragma omp section
    protW00tW->evaluate(MMt);
#pragma omp section
    prot00WW0->evaluate(MMt);
#pragma omp section
    prot000  ->evaluate(MMt);
#pragma omp section
    prot0W00 ->evaluate(MMt);

#pragma omp section
    protH0tt0->evaluate(MMt);
#pragma omp section
    protH0t00->evaluate(MMt);
#pragma omp section
    prot0Htt0->evaluate(MMt);
#pragma omp section
    prot0H0t0->evaluate(MMt);
#pragma omp section
    prot00ttH->evaluate(MMt);
#pragma omp section
    protHtt0t->evaluate(MMt);
#pragma omp section
    prot00t00->evaluate(MMt);
#pragma omp section
    prot000t0->evaluate(MMt);
    }
#else
    protWt000->evaluate(MMt);
    prot0ttHt->evaluate(MMt);
    prot0ttZt->evaluate(MMt);
    prot0tt0t->evaluate(MMt);
    prottH0H->evaluate(MMt);
    prottZ0Z->evaluate(MMt);
    
    protHHttH->evaluate(MMt);
    protHZttZ->evaluate(MMt);
    protHWt0W->evaluate(MMt);
    protHttHt->evaluate(MMt);
    protHttZt->evaluate(MMt);
    protZZttH->evaluate(MMt);
    protZWt0W->evaluate(MMt);
    protZttZt->evaluate(MMt);
    protZ0tW0->evaluate(MMt);
    protWW00Z->evaluate(MMt);
    protW00tW->evaluate(MMt);
    prot00WW0->evaluate(MMt);
    prot000  ->evaluate(MMt);
    prot0W00 ->evaluate(MMt);
    
    protH0tt0->evaluate(MMt);
    protH0t00->evaluate(MMt);
    prot0Htt0->evaluate(MMt);
    prot0H0t0->evaluate(MMt);
    prot00ttH->evaluate(MMt);
    protHtt0t->evaluate(MMt);
    prot00t00->evaluate(MMt);
    prot000t0->evaluate(MMt);
#endif
    t.elapsed();
  }


  tt<MS>::tt(MSinput sm, double mu2_)
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
  
    protWt000 = std::unique_ptr<Tsil>(new Tsil(mmW, mmt,   0,   0,   0, mu2) );
    prot0ttHt = std::unique_ptr<Tsil>(new Tsil(0  , mmt, mmt, mmH, mmt, mu2) );
    prot0ttZt = std::unique_ptr<Tsil>(new Tsil(0  , mmt, mmt, mmZ, mmt, mu2) );
    prot0tt0t = std::unique_ptr<Tsil>(new Tsil(0  , mmt, mmt,   0, mmt, mu2) );
    prottH0H =  std::unique_ptr<TsilSTU>(new TsilSTU ( mmt, mmH,   0, mmH, mu2) );
    prottZ0Z =  std::unique_ptr<TsilSTU>(new TsilSTU ( mmt, mmZ,   0, mmZ, mu2) );
  
    protHHttH = std::unique_ptr<Tsil>(new Tsil(mmH, mmH, mmt, mmt, mmH, mu2) );
    protHZttZ = std::unique_ptr<Tsil>(new Tsil(mmH, mmZ, mmt, mmt, mmZ, mu2) );
    protHWt0W = std::unique_ptr<Tsil>(new Tsil(mmH, mmW, mmt,   0, mmW, mu2) );
    protHttHt = std::unique_ptr<Tsil>(new Tsil(mmH, mmt, mmt, mmH, mmt, mu2) );
    protHttZt = std::unique_ptr<Tsil>(new Tsil(mmH, mmt, mmt, mmZ, mmt, mu2) );
    protZZttH = std::unique_ptr<Tsil>(new Tsil(mmZ, mmZ, mmt, mmt, mmH, mu2) );
    protZWt0W = std::unique_ptr<Tsil>(new Tsil(mmZ, mmW, mmt,   0, mmW, mu2) );
    protZttZt = std::unique_ptr<Tsil>(new Tsil(mmZ, mmt, mmt, mmZ, mmt, mu2) );
    protZ0tW0 = std::unique_ptr<Tsil>(new Tsil(mmZ,   0, mmt, mmW,   0, mu2) );
    protWW00Z = std::unique_ptr<Tsil>(new Tsil(mmW, mmW,   0,   0, mmZ, mu2) );
    protW00tW = std::unique_ptr<Tsil>(new Tsil(mmW,   0,   0, mmt, mmW, mu2) );
    prot00WW0 = std::unique_ptr<Tsil>(new Tsil(0  ,   0, mmW, mmW,   0, mu2) );
    prot000   = std::unique_ptr<TsilST>(new TsilST(0,             0,   0, mu2) );
    prot0W00  = std::unique_ptr<TsilSTU>(new TsilSTU(0,     mmW,   0,   0, mu2) );

    protH0tt0 = std::unique_ptr<Tsil>(new Tsil(mmH,   0, mmt, mmt,   0, mu2) );
    protH0t00 = std::unique_ptr<Tsil>(new Tsil(mmH,   0, mmt,   0,   0, mu2) );
    prot0Htt0 = std::unique_ptr<Tsil>(new Tsil(  0, mmH, mmt, mmt,   0, mu2) );
    prot0H0t0 = std::unique_ptr<Tsil>(new Tsil(  0, mmH,   0, mmt,   0, mu2) );
    prot00ttH = std::unique_ptr<Tsil>(new Tsil(  0,   0, mmt, mmt, mmH, mu2) );
    protHtt0t = std::unique_ptr<Tsil>(new Tsil(mmH, mmt, mmt,   0, mmt, mu2) );
    prot00t00 = std::unique_ptr<Tsil>(new Tsil(  0,   0, mmt,   0,   0, mu2) );
    prot000t0 = std::unique_ptr<Tsil>(new Tsil(  0,   0,   0, mmt,   0, mu2) );

  
    Timer t;
#ifdef _OPENMP
#pragma omp parallel sections
    {
#pragma omp section
      protWt000->evaluate(mmt);
#pragma omp section
      prot0ttHt->evaluate(mmt);
#pragma omp section
      prot0ttZt->evaluate(mmt);
#pragma omp section
      prot0tt0t->evaluate(mmt);
#pragma omp section
      prottH0H->evaluate(mmt);
#pragma omp section
      prottZ0Z->evaluate(mmt);

#pragma omp section
      protHHttH->evaluate(mmt);
#pragma omp section
      protHZttZ->evaluate(mmt);
#pragma omp section
      protHWt0W->evaluate(mmt);
#pragma omp section
      protHttHt->evaluate(mmt);
#pragma omp section
      protHttZt->evaluate(mmt);
#pragma omp section
      protZZttH->evaluate(mmt);
#pragma omp section
      protZWt0W->evaluate(mmt);
#pragma omp section
      protZttZt->evaluate(mmt);
#pragma omp section
      protZ0tW0->evaluate(mmt);
#pragma omp section
      protWW00Z->evaluate(mmt);
#pragma omp section
      protW00tW->evaluate(mmt);
#pragma omp section
      prot00WW0->evaluate(mmt);
#pragma omp section
      prot000  ->evaluate(mmt);
#pragma omp section
      prot0W00 ->evaluate(mmt);

#pragma omp section
      protH0tt0->evaluate(mmt);
#pragma omp section
      protH0t00->evaluate(mmt);
#pragma omp section
      prot0Htt0->evaluate(mmt);
#pragma omp section
      prot0H0t0->evaluate(mmt);
#pragma omp section
      prot00ttH->evaluate(mmt);
#pragma omp section
      protHtt0t->evaluate(mmt);
#pragma omp section
      prot00t00->evaluate(mmt);
#pragma omp section
      prot000t0->evaluate(mmt);
    }
#else
      protWt000->evaluate(mmt);
      prot0ttHt->evaluate(mmt);
      prot0ttZt->evaluate(mmt);
      prot0tt0t->evaluate(mmt);
      prottH0H->evaluate(mmt);
      prottZ0Z->evaluate(mmt);
      
      protHHttH->evaluate(mmt);
      protHZttZ->evaluate(mmt);
      protHWt0W->evaluate(mmt);
      protHttHt->evaluate(mmt);
      protHttZt->evaluate(mmt);
      protZZttH->evaluate(mmt);
      protZWt0W->evaluate(mmt);
      protZttZt->evaluate(mmt);
      protZ0tW0->evaluate(mmt);
      protWW00Z->evaluate(mmt);
      protW00tW->evaluate(mmt);
      prot00WW0->evaluate(mmt);
      prot000  ->evaluate(mmt);
      prot0W00 ->evaluate(mmt);
      
      protH0tt0->evaluate(mmt);
      protH0t00->evaluate(mmt);
      prot0Htt0->evaluate(mmt);
      prot0H0t0->evaluate(mmt);
      prot00ttH->evaluate(mmt);
      protHtt0t->evaluate(mmt);
      prot00t00->evaluate(mmt);
      prot000t0->evaluate(mmt);
#endif
      t.elapsed();
  
  }
} // namespace mr
