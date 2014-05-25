// #include <omp.h>
#include <tt.hpp>
#include "timer.hpp"


tt::tt(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
  MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
{
  init(MMW, MMZ, MMH, MMt, mu2);
}

tt::tt(SMinput sm, long double mu2_)
{
  MMW = sm.MMW();
  MMZ = sm.MMZ();
  MMH = sm.MMH();
  MMt = sm.MMt();
  mu2 = mu2_;

  init(sm.MMW(), sm.MMZ(), sm.MMH(), sm.MMt(), mu2_);
}


void tt::init(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_)
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
    protos[i]->evaluate(MMt);
  t2.elapsed();

}

const long double tt::EPAIR2;
