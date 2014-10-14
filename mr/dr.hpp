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

// alpha*alphaS
// Tsil* protWt000;
// Tsil* prot0ttHt;
// Tsil* prot0ttZt;
// Tsil* prot0tt0t;
// TsilSTU* prottH0H;
// TsilSTU* prottZ0Z;
//   // alpha^2
//   Tsil* protHHttH;
//   Tsil* protHZttZ;
//   Tsil* protHWt0W;
//   Tsil* protHttHt;
//   Tsil* protHttZt;
//   Tsil* protZZttH;
//   Tsil* protZWt0W;
//   Tsil* protZttZt;
//   Tsil* protZ0tW0;
//   Tsil* protWW00Z;
//   Tsil* protW00tW;
//   Tsil* prot00WW0;
//   TsilSTU* prot0W00;
//   TsilST* prot000;

//   // Gaugeless limit
//   Tsil* protH0tt0;
//   Tsil* protH0t00;
//   Tsil* prot0Htt0;
//   Tsil* prot0H0t0;
//   Tsil* prot00ttH;
//   Tsil* protHtt0t;
//   Tsil* prot00t00;
//   Tsil* prot000t0;

// TsilST* protos[28];

static const long double EPAIR2 = -1.; 

public:
dr()
{
}

  dr(long double,long double,long double,long double,long double);
  
  dr(SMinput, long double);
  
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
