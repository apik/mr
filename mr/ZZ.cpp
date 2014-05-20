// #include <omp.h>
#include <ZZ.hpp>
#include "timer.hpp"

ZZ::ZZ(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_):
  MMW(MMW_), MMZ(MMZ_), MMH(MMH_), MMt(MMt_), mu2(mu2_)
{
  init(MMW, MMZ, MMH, MMt, mu2);
}

ZZ::ZZ(SMinput sm, long double mu2_)
{
  MMW = sm.MMW();
  MMZ = sm.MMZ();
  MMH = sm.MMH();
  MMt = sm.MMt();
  mu2 = mu2_;

  init(sm.MMW(), sm.MMZ(), sm.MMH(), sm.MMt(), mu2_);
}


void ZZ::init(long double MMW_,long double MMZ_,long double MMH_,long double MMt_,long double mu2_)
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
