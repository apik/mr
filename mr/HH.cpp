// #include <omp.h>
#include <HH.hpp>
#include "timer.hpp"


// HH::HH(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
//   MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
// {
//   init(MMW, MMZ, MMH, MMt, mu2);
// }

HH::HH(SMinput sm, long double mu2_)
{
  MMb = sm.MMb();
  MMW = sm.MMW();
  MMZ = sm.MMZ();
  MMH = sm.MMH();
  MMt = sm.MMt();
  mu2 = mu2_;

  init();
}


void HH::init()
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
