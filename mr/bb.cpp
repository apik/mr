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

#include <bb.hpp>
#include "timer.hpp"


namespace mr
{
  bb<OS>::bb(long double MMb_, long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
    MMb(MMb_), MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
  {
    init(MMb, MMW, MMZ, MMH, MMt, mu2);
  }

  bb<OS>::bb(OSinput sm, long double mu2_)
  {
    MMb = sm.MMb();
    MMW = sm.MMW();
    MMZ = sm.MMZ();
    MMH = sm.MMH();
    MMt = sm.MMt();
    mu2 = mu2_;

    init(sm.MMb(),sm.MMW(), sm.MMZ(), sm.MMH(), sm.MMt(), mu2_);
  }


  void bb<OS>::init(long double MMb_, long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_)
  {
  
    CW = sqrt(MMW/MMZ);
    SW = sqrt(1-MMW/MMZ);
  
    prot0bb0b = new Tsil(   0, MMb, MMb,   0, MMb, mu2);

    // Timer t2;
    prot0bb0b->evaluate(MMb);
    // t2.elapsed();

  }

  bb<MS>::bb(MSinput sm, long double mu2_)
  {
    mmb = sm.mmb();
    mmW = sm.mmW();
    mmZ = sm.mmZ();
    mmH = sm.mmH();
    mmt = sm.mmt();
    mu2 = mu2_;

    init();
  }


  void bb<MS>::init()
  {
  
    c = sqrt(mmW/mmZ);
    s = sqrt(1-mmW/mmZ);

    prot0bb0b = new Tsil(   0, mmb, mmb,   0, mmb, mu2);
    
    // Timer t2;
    prot0bb0b->evaluate(mmb);
    // t2.elapsed();
  }
} // namespace mr
