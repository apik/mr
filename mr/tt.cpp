#include <omp.h>
#include <tt.hpp>
#include "timer.hpp"
tt::tt(long double MMt_,long double MMH_,long double MMW_,long double MMZ_,long double mmu_):
  MMt(MMt_), MMH(MMH_), MMW(MMW_), MMZ(MMZ_), mmu(mmu_)
{
  // protWt000 = new Tsil(MMW, MMt,   0,   0,   0, mmu, MMt);
  // prot0ttHt = new Tsil(0  , MMt, MMt, MMH, MMt, mmu, MMt);
  // prot0ttZt = new Tsil(0  , MMt, MMt, MMZ, MMt, mmu, MMt);
  // prot0tt0t = new Tsil(0  , MMt, MMt,   0, MMt, mmu, MMt);
  // prottH0H =  new TsilSTU ( MMt, MMH,   0, MMH, mmu, MMt);
  // prottZ0Z =  new TsilSTU ( MMt, MMZ,   0, MMZ, mmu, MMt);
  
  // protHHttH = new Tsil(MMH, MMH, MMt, MMt, MMH, mmu, MMt);
  // protHZttZ = new Tsil(MMH, MMZ, MMt, MMt, MMZ, mmu, MMt);
  // protHWt0W = new Tsil(MMH, MMW, MMt,   0, MMW, mmu, MMt);
  // protHttHt = new Tsil(MMH, MMt, MMt, MMH, MMt, mmu, MMt);
  // protHttZt = new Tsil(MMH, MMt, MMt, MMZ, MMt, mmu, MMt);
  // protZZttH = new Tsil(MMZ, MMZ, MMt, MMt, MMH, mmu, MMt);
  // protZWt0W = new Tsil(MMZ, MMW, MMt,   0, MMW, mmu, MMt);
  // protZttZt = new Tsil(MMZ, MMt, MMt, MMZ, MMt, mmu, MMt);
  // protZ0tW0 = new Tsil(MMZ,   0, MMt, MMW,   0, mmu, MMt);
  // protWW00Z = new Tsil(MMW, MMW,   0,   0, MMZ, mmu, MMt);
  // protW00tW = new Tsil(MMW,   0,   0, MMt, MMW, mmu, MMt);
  // prot00WW0 = new Tsil(0  ,   0, MMW, MMW,   0, mmu, MMt);
  // prot000   = new TsilST(0  , 0 ,0 , mmu, MMt);
  // prot0W00  = new TsilSTU(0,MMW,0,0,mmu,MMt);


  protos[0]  = protWt000 = new Tsil(MMW, MMt,   0,   0,   0, mmu);
  protos[1]  = prot0ttHt = new Tsil(0  , MMt, MMt, MMH, MMt, mmu);
  protos[2]  = prot0ttZt = new Tsil(0  , MMt, MMt, MMZ, MMt, mmu);
  protos[3]  = prot0tt0t = new Tsil(0  , MMt, MMt,   0, MMt, mmu);
  protos[4]  = prottH0H =  new TsilSTU ( MMt, MMH,   0, MMH, mmu);
  protos[5]  = prottZ0Z =  new TsilSTU ( MMt, MMZ,   0, MMZ, mmu);
  
  protos[6]  = protHHttH = new Tsil(MMH, MMH, MMt, MMt, MMH, mmu);
  protos[7]  = protHZttZ = new Tsil(MMH, MMZ, MMt, MMt, MMZ, mmu);
  protos[8]  = protHWt0W = new Tsil(MMH, MMW, MMt,   0, MMW, mmu);
  protos[9]  = protHttHt = new Tsil(MMH, MMt, MMt, MMH, MMt, mmu);
  protos[10] = protHttZt = new Tsil(MMH, MMt, MMt, MMZ, MMt, mmu);
  protos[11] = protZZttH = new Tsil(MMZ, MMZ, MMt, MMt, MMH, mmu);
  protos[12] = protZWt0W = new Tsil(MMZ, MMW, MMt,   0, MMW, mmu);
  protos[13] = protZttZt = new Tsil(MMZ, MMt, MMt, MMZ, MMt, mmu);
  protos[14] = protZ0tW0 = new Tsil(MMZ,   0, MMt, MMW,   0, mmu);
  protos[15] = protWW00Z = new Tsil(MMW, MMW,   0,   0, MMZ, mmu);
  protos[16] = protW00tW = new Tsil(MMW,   0,   0, MMt, MMW, mmu);
  protos[17] = prot00WW0 = new Tsil(0  ,   0, MMW, MMW,   0, mmu);
  protos[18] = prot000   = new TsilST(0,             0,   0, mmu);
  protos[19] = prot0W00  = new TsilSTU(0,     MMW,   0,   0, mmu);

  Timer t1;

  int TID = 0;
  omp_set_num_threads(10);
#pragma omp parallel private(TID)
  {
    TID = omp_get_thread_num();
    std::cout << "Evaluating proto [" << TID << "]" <<  std::endl;
    protos[TID]->evaluate(MMt);
    
  }
  
  t1.elapsed();

  Timer t2;
  for(int i = 0 ; i < 20; i++)
    protos[i]->evaluate(MMt);
  t2.elapsed();

  std::cout << "Constr~!!\n";
}
