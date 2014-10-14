// #include <omp.h>
#include <WW.hpp>
#include "timer.hpp"

// WW::WW(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
//   MMb(MMb_), MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
// {
//   init(MMB, MMW, MMZ, MMH, MMt, mu2);
// }

WW::WW(SMinput sm, long double mu2_)
{
  MMb = sm.MMb();
  MMW = sm.MMW();
  MMZ = sm.MMZ();
  MMH = sm.MMH();
  MMt = sm.MMt();
  mu2 = mu2_;

  init();
}


void WW::init()
{

  CW = sqrt(MMW/MMZ);
  SW = sqrt(1-MMW/MMZ);

  protos[0] =  protWHHWW = new Tsil(MMW, MMH, MMH, MMW, MMW, mu2);
  protos[1] =  protWHZWW = new Tsil(MMW, MMH, MMZ, MMW, MMW, mu2);
  protos[2] =  protWZZWW = new Tsil(MMW, MMZ, MMZ, MMW, MMW, mu2);
  protos[3] =  protWWHHH = new Tsil(MMW, MMW, MMH, MMH, MMH, mu2);
  protos[4] =  protWWHZZ = new Tsil(MMW, MMW, MMH, MMZ, MMZ, mu2);
  protos[5] =  protWWZZH = new Tsil(MMW, MMW, MMZ, MMZ, MMH, mu2);
  protos[6] =  protWtZ00 = new Tsil(MMW, MMt, MMZ,   0,   0, mu2);
  protos[7] =  protW0HWW = new Tsil(MMW,   0, MMH, MMW, MMW, mu2);
  protos[8] =  protW0Htt = new Tsil(MMW,   0, MMH, MMt, MMt, mu2);
  protos[9] =  protW0ZWW = new Tsil(MMW,   0, MMZ, MMW, MMW, mu2);
  protos[10] =  protW0Ztt = new Tsil(MMW,   0, MMZ, MMt, MMt, mu2);
  protos[11] =  protW0Z00 = new Tsil(MMW,   0, MMZ,   0,   0, mu2);
  protos[12] =  prot0WW0W = new Tsil(  0, MMW, MMW,   0, MMW, mu2);
  protos[13] =  prot0Wt0t = new Tsil(  0, MMW, MMt,   0, MMt, mu2);
  protos[14] =  prot0W0Z0 = new Tsil(  0, MMW,   0, MMZ,   0, mu2);
  protos[15] =  prot00Wt0 = new Tsil(  0,   0, MMW, MMt,   0, mu2);
  protos[16] =  prot00W00 = new Tsil(  0,   0, MMW,   0,   0, mu2);
  protos[17] =  prot00ttZ = new Tsil(  0,   0, MMt, MMt, MMZ, mu2);
  protos[18] =  prot00tt0 = new Tsil(  0,   0, MMt, MMt,   0, mu2);
  protos[19] =  prot0000Z = new Tsil(  0,   0,   0,  0, MMZ, mu2);
  protos[20] =  prot00000 = new Tsil(  0,   0,   0,   0,  0, mu2);
  protos[21] =  protWH0H = new TsilSTU(MMW, MMH,   0, MMH, mu2);
  protos[22] =  protWZ0Z = new TsilSTU(MMW, MMZ,   0, MMZ, mu2);
  protos[23] =  protHW00 = new TsilSTU(MMH, MMW,   0,   0, mu2);



  // protos[19] = prot0W00  = new TsilSTU(0,     MMW,   0,   0, mu2);

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
  for(int i = 0 ; i < 24; i++)
    protos[i]->evaluate(MMW);
  t2.elapsed();

}

const long double WW::EPAIR2;
