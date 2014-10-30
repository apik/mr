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
#include <bb.hpp>
#include "timer.hpp"


bb::bb(long double MMb_, long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
  MMb(MMb_), MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
{
  init(MMb, MMW, MMZ, MMH, MMt, mu2);
}

bb::bb(OSinput sm, long double mu2_)
{
  MMb = sm.MMb();
  MMW = sm.MMW();
  MMZ = sm.MMZ();
  MMH = sm.MMH();
  MMt = sm.MMt();
  mu2 = mu2_;

  init(sm.MMb(),sm.MMW(), sm.MMZ(), sm.MMH(), sm.MMt(), mu2_);
}


void bb::init(long double MMb_, long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_)
{
  
  CW = sqrt(MMW/MMZ);
  SW = sqrt(1-MMW/MMZ);
  
  prot0bb0b = new Tsil(   0, MMb, MMb,   0, MMb, mu2);
//   Timer t1;

//   int TID = 0;
//   omp_set_num_threads(10);
// #pragma omp parallel private(TID)
//   {
//     TID = omp_get_thread_num();
//     std::cout << "Evaluating proto [" << TID << "]" <<  std::endl;
//     protos[TID]->evaluate(MMt);
    
//   }
  
//   t1.elapsed();

  Timer t2;
  // for(int i = 0 ; i < 20; i++)
  //   protos[i]->evaluate(MMt);

  prot0bb0b->evaluate(MMb);
  t2.elapsed();

}

const long double bb::EPAIR2;
