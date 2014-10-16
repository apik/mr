#ifndef __DR_HPP__
#define __DR_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"


class dr
{

long double MMb,MMt, MMH, MMW, MMZ, mu2;
long double SW, CW;


public:
dr()
{
}

  dr(long double,long double,long double,long double,long double);
  
  dr(OSinput, long double);
  
  void init();


  
  
  

  // 
  // Pure QCD part m_ij=mY_ij by definition
  // 
  // std::complex<long double> m01();

  // std::complex<long double> m02(size_t nL = 5);

  // std::complex<long double> m03(size_t nL = 5);

  // 
  // Mass corrections
  // 
  // std::complex<long double> m10(size_t nL = 2, size_t nH = 1);
  
  // std::complex<long double> m11(size_t nL = 2, size_t nH = 1);

  // std::complex<long double> m20(size_t nL = 2, size_t nH = 1);


  // 
  // Delta-r
  // 
  std::complex<long double> dr10(size_t nL = 2, size_t nH = 1);
  
  std::complex<long double> dr11(size_t nL = 2, size_t nH = 1);

  std::complex<long double> dr20(size_t nL = 2, size_t nH = 1);


  // Gaugeless limit
  std::complex<long double> drgl10(size_t nL = 2, size_t nH = 1);
  
  std::complex<long double> drgl11(size_t nL = 2, size_t nH = 1);

  std::complex<long double> drgl20(size_t nL = 2, size_t nH = 1);

  
};


#endif  //  __DR_HPP__
