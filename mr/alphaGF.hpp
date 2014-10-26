#ifndef __ALPHAGF_HPP__
#define __ALPHAGF_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"
#include "base.hpp"

class alphaGF
{

  long double MMb, MMt, MMH, MMW, MMZ, mu2;
  long double SW, CW;

  Tsil* protWHHWW;
  Tsil* protWHZWW;
  Tsil* protWZZWW;
  Tsil* protWWHHH;
  Tsil* protWWHZZ;
  Tsil* protWWZZH;
  Tsil* protWtZ00;
  Tsil* protW0HWW;
  Tsil* protW0Htt;
  Tsil* protW0ZWW;
  Tsil* protW0Ztt;
  Tsil* protW0Z00;
  Tsil* prot0WW0W;
  Tsil* prot0Wt0t;
  Tsil* prot0W0Z0;
  Tsil* prot00Wt0;
  Tsil* prot00W00;
  Tsil* prot00ttZ;
  Tsil* prot00tt0;
  Tsil* prot0000Z;
  Tsil* prot00000;
  // ZZ part

  Tsil* protZHHZZ;
  Tsil* protZZHHH;
  Tsil* protZWHWW;
  Tsil* prottZtHt;
  Tsil* protWWWWH;
  Tsil* protWWWWZ;
  Tsil* protWWWW0;
  Tsil* protWtWt0;
  Tsil* protW0W0t;
  Tsil* protW0W00;
  Tsil* protttttH;
  Tsil* protttttZ;
  Tsil* prottttt0;
  Tsil* prott0t0W;
  // Tsil* prot0000Z;
  Tsil* prot0000W;
  // Tsil* prot00000;
  Tsil* protWZWHW;
  TsilSTU* protHZ00;


  TsilSTU* protWH0H;
  TsilSTU* protWZ0Z;
  TsilSTU* protHW00;
  
  TsilST* protos[24];
public:
  alphaGF()
  {
  }

  alphaGF(long double,long double,long double,long double,long double);

  alphaGF(OSinput, long double);

  void init();
  
 
  std::complex<long double> a10(size_t nL = 2, size_t nH = 1);

  
  std::complex<long double> a11(size_t nL = 2, size_t nH = 1);
  

  std::complex<long double> a20(size_t nL = 2, size_t nH = 1);
};


#endif  // __ALPHAGF_HPP__
