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
#include <HH.hpp>
#include "timer.hpp"


// HH::HH(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
//   MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
// {
//   init(MMW, MMZ, MMH, MMt, mu2);
// }

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
  
  protos[0] = protHHHHH = new Tsil(MMH, MMH, MMH, MMH, MMH, mu2);
  protos[1] = protHZHZZ = new Tsil(MMH, MMZ, MMH, MMZ, MMZ, mu2);
  protos[2] = protHWHWW = new Tsil(MMH, MMW, MMH, MMW, MMW, mu2);
  protos[3] = protHtHtt = new Tsil(MMH, MMt, MMH, MMt, MMt, mu2);
  protos[4] = protZZZZH = new Tsil(MMZ, MMZ, MMZ, MMZ, MMH, mu2);
  protos[5] = protZWZWW = new Tsil(MMZ, MMW, MMZ, MMW, MMW, mu2);
  protos[6] = protZtZtt = new Tsil(MMZ, MMt, MMZ, MMt, MMt, mu2);
  protos[7] = protWWWWH = new Tsil(MMW, MMW, MMW, MMW, MMH, mu2);
  protos[8] = protWWWWZ = new Tsil(MMW, MMW, MMW, MMW, MMZ, mu2);
  protos[9] = protWWWW0 = new Tsil(MMW, MMW, MMW, MMW,   0, mu2);
  protos[10] = protWtWt0 = new Tsil(MMW, MMt, MMW, MMt,   0, mu2);
  protos[11] = protttttH = new Tsil(MMt, MMt, MMt, MMt, MMH, mu2);
  protos[12] = protttttZ = new Tsil(MMt, MMt, MMt, MMt, MMZ, mu2);
  protos[13] = prottttt0 = new Tsil(MMt, MMt, MMt, MMt,   0, mu2);
  protos[14] = protZZ00 = new TsilSTU(MMZ, MMZ,    0,   0, mu2);
  protos[15] = protWW00 = new TsilSTU(MMW, MMW,    0,   0, mu2);

  protos[16] = prot0H0H0 = new Tsil(  0, MMH,   0, MMH,   0, mu2);
  protos[17] = prot0t0tt = new Tsil(  0, MMt,   0, MMt, MMt, mu2);
  protos[18] = prot0t0t0 = new Tsil(  0, MMt,   0, MMt,   0, mu2);
  protos[19] = prot0000H = new Tsil(  0,   0,   0,   0, MMH, mu2);




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
  for(int i = 0 ; i < 20; i++)
    protos[i]->evaluate(MMH);
  t2.elapsed();

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
  
  protos[0] = protHHHHH = new Tsil(mmH, mmH, mmH, mmH, mmH, mu2);
  protos[1] = protHZHZZ = new Tsil(mmH, mmZ, mmH, mmZ, mmZ, mu2);
  protos[2] = protHWHWW = new Tsil(mmH, mmW, mmH, mmW, mmW, mu2);
  protos[3] = protHtHtt = new Tsil(mmH, mmt, mmH, mmt, mmt, mu2);
  protos[4] = protZZZZH = new Tsil(mmZ, mmZ, mmZ, mmZ, mmH, mu2);
  protos[5] = protZWZWW = new Tsil(mmZ, mmW, mmZ, mmW, mmW, mu2);
  protos[6] = protZtZtt = new Tsil(mmZ, mmt, mmZ, mmt, mmt, mu2);
  protos[7] = protWWWWH = new Tsil(mmW, mmW, mmW, mmW, mmH, mu2);
  protos[8] = protWWWWZ = new Tsil(mmW, mmW, mmW, mmW, mmZ, mu2);
  protos[9] = protWWWW0 = new Tsil(mmW, mmW, mmW, mmW,   0, mu2);
  protos[10] = protWtWt0 = new Tsil(mmW, mmt, mmW, mmt,   0, mu2);
  protos[11] = protttttH = new Tsil(mmt, mmt, mmt, mmt, mmH, mu2);
  protos[12] = protttttZ = new Tsil(mmt, mmt, mmt, mmt, mmZ, mu2);
  protos[13] = prottttt0 = new Tsil(mmt, mmt, mmt, mmt,   0, mu2);
  protos[14] = protZZ00 = new TsilSTU(mmZ, mmZ,    0,   0, mu2);
  protos[15] = protWW00 = new TsilSTU(mmW, mmW,    0,   0, mu2);

  protos[16] = prot0H0H0 = new Tsil(  0, mmH,   0, mmH,   0, mu2);
  protos[17] = prot0t0tt = new Tsil(  0, mmt,   0, mmt, mmt, mu2);
  protos[18] = prot0t0t0 = new Tsil(  0, mmt,   0, mmt,   0, mu2);
  protos[19] = prot0000H = new Tsil(  0,   0,   0,   0, mmH, mu2);




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
  for(int i = 0 ; i < 20; i++)
    protos[i]->evaluate(mmH);
  t2.elapsed();

}
