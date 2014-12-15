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

  protos[20]  = protH0tt0 = new Tsil(MMH,   0, MMt, MMt,   0, mu2);
  protos[21]  = protH0t00 = new Tsil(MMH,   0, MMt,   0,   0, mu2);
  protos[22]  = prot0Htt0 = new Tsil(  0, MMH, MMt, MMt,   0, mu2);
  protos[23]  = prot0H0t0 = new Tsil(  0, MMH,   0, MMt,   0, mu2);
  protos[24]  = prot00ttH = new Tsil(  0,   0, MMt, MMt, MMH, mu2);
  protos[25]  = protHtt0t = new Tsil(MMH, MMt, MMt,   0, MMt, mu2);
  protos[26]  = prot00t00 = new Tsil(  0,   0, MMt,   0,   0, mu2);
  protos[27]  = prot000t0 = new Tsil(  0,   0,   0, MMt,   0, mu2);

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
  for(int i = 0 ; i < 28; i++)
    protos[i]->evaluate(MMt);
  t2.elapsed();

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
  for(int i = 0 ; i < 28; i++)
    protos[i]->evaluate(mmt);
  t2.elapsed();

}



std::complex<long double> DIA(int i)
{
  // return i==1?1:0;
  return 0;
}
 
std::pair<long double,long double> tt<OS>::test2(long double epsabs,long double epsrel)
{

  std::cout << "Test 2" << std::endl;

  long double DH = 1-MMH/MMt;
  std::vector<std::complex<long double> > rexact(139);
  std::vector<std::complex<long double> > rexpan(139);

// #include "dump.hpp"

  long double accEX = 0;
  long double accEP = 0;

  for(int i = 0; i < 139; i++)
    {
      const size_t fw = 15;

      long double diff  = fabs((rexact[i] - rexpan[i]).real());
      long double exact = fabs(rexact[i].real());

      // if (true ||diff > epsabs || diff/exact > epsrel)
      //   {
      std::cout << std::scientific << std::setprecision(fw-10);
      std::cout << "Dia(" << i << ") ";
      std::cout << "   " << std::setw(fw) <<  rexact[i].real()
                << "   " << std::setw(fw) <<  rexpan[i].real()
                << "   " << std::setw(fw) << diff <<std::endl; 
      // }

      accEX += rexact[i].real();
      accEP += rexpan[i].real();

    }
  return std::make_pair(accEX,accEP);
}



// const long double tt::EPAIR2;

