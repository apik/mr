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
#include <ZZ.hpp>
#include "timer.hpp"

// ZZ::ZZ(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
//   MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
// {
//   init(MMW, MMZ, MMH, MMt, mu2);
// }

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

  protos[0] = protZHHZZ = new Tsil(MMZ, MMH, MMH, MMZ, MMZ, mu2);
  protos[1] = protZZHHH = new Tsil(MMZ, MMZ, MMH, MMH, MMH, mu2);
  protos[2] = protZWHWW = new Tsil(MMZ, MMW, MMH, MMW, MMW, mu2);
  protos[3] = prottZtHt = new Tsil(MMt, MMZ, MMt, MMH, MMt, mu2);
  protos[4] = protWWWWH = new Tsil(MMW, MMW, MMW, MMW, MMH, mu2);
  protos[5] = protWWWWZ = new Tsil(MMW, MMW, MMW, MMW, MMZ, mu2);
  protos[6] = protWWWW0 = new Tsil(MMW, MMW, MMW, MMW,   0, mu2);
  protos[7] = protWtWt0 = new Tsil(MMW, MMt, MMW, MMt,   0, mu2);
  protos[8] = protW0W0t = new Tsil(MMW,   0, MMW,   0, MMt, mu2);
  protos[9] = protW0W00 = new Tsil(MMW,   0, MMW,   0,   0, mu2);
  protos[10] = protttttH = new Tsil(MMt, MMt, MMt, MMt, MMH, mu2);
  protos[11] = protttttZ = new Tsil(MMt, MMt, MMt, MMt, MMZ, mu2);
  protos[12] = prottttt0 = new Tsil(MMt, MMt, MMt, MMt,   0, mu2);
  protos[13] = prott0t0W = new Tsil(MMt,   0, MMt,   0, MMW, mu2);
  protos[14] = prot0000Z = new Tsil(  0,   0,   0,   0, MMZ, mu2);
  protos[15] = prot0000W = new Tsil(  0,   0,   0,   0, MMW, mu2);
  protos[16] = prot00000 = new Tsil(  0,   0,   0,   0,   0, mu2);
  protos[17] = protWZWHW = new Tsil(MMW, MMZ, MMW, MMH, MMW, mu2);
  protos[18] = protHZ00  = new TsilSTU(MMH, MMZ,   0,  0, mu2);


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
  for(int i = 0 ; i < 19; i++)
    protos[i]->evaluate(MMZ);
  t2.elapsed();

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

  protos[0] = protZHHZZ = new Tsil(mmZ, mmH, mmH, mmZ, mmZ, mu2);
  protos[1] = protZZHHH = new Tsil(mmZ, mmZ, mmH, mmH, mmH, mu2);
  protos[2] = protZWHWW = new Tsil(mmZ, mmW, mmH, mmW, mmW, mu2);
  protos[3] = prottZtHt = new Tsil(mmt, mmZ, mmt, mmH, mmt, mu2);
  protos[4] = protWWWWH = new Tsil(mmW, mmW, mmW, mmW, mmH, mu2);
  protos[5] = protWWWWZ = new Tsil(mmW, mmW, mmW, mmW, mmZ, mu2);
  protos[6] = protWWWW0 = new Tsil(mmW, mmW, mmW, mmW,   0, mu2);
  protos[7] = protWtWt0 = new Tsil(mmW, mmt, mmW, mmt,   0, mu2);
  protos[8] = protW0W0t = new Tsil(mmW,   0, mmW,   0, mmt, mu2);
  protos[9] = protW0W00 = new Tsil(mmW,   0, mmW,   0,   0, mu2);
  protos[10] = protttttH = new Tsil(mmt, mmt, mmt, mmt, mmH, mu2);
  protos[11] = protttttZ = new Tsil(mmt, mmt, mmt, mmt, mmZ, mu2);
  protos[12] = prottttt0 = new Tsil(mmt, mmt, mmt, mmt,   0, mu2);
  protos[13] = prott0t0W = new Tsil(mmt,   0, mmt,   0, mmW, mu2);
  protos[14] = prot0000Z = new Tsil(  0,   0,   0,   0, mmZ, mu2);
  protos[15] = prot0000W = new Tsil(  0,   0,   0,   0, mmW, mu2);
  protos[16] = prot00000 = new Tsil(  0,   0,   0,   0,   0, mu2);
  protos[17] = protWZWHW = new Tsil(mmW, mmZ, mmW, mmH, mmW, mu2);
  protos[18] = protHZ00  = new TsilSTU(mmH, mmZ,   0,  0, mu2);


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
  for(int i = 0 ; i < 19; i++)
    protos[i]->evaluate(mmZ);
  t2.elapsed();

}

