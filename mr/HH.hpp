#ifndef __HH_HPP__
#define __HH_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"


class HH
{

  long double MMt, MMH, MMW, MMZ, mu2;
  long double SW, CW;
  
  Tsil* protHHHHH;
  Tsil* protHZHZZ;
  Tsil* protHWHWW;
  Tsil* protHtHtt;
  Tsil* protZZZZH;
  Tsil* protZWZWW;
  Tsil* protZtZtt;
  Tsil* protWWWWH;
  Tsil* protWWWWZ;
  Tsil* protWWWW0;
  Tsil* protWtWt0;
  Tsil* protttttH;
  Tsil* protttttZ;
  Tsil* prottttt0;
  TsilSTU* protZZ00;
  TsilSTU* protWW00;
  
  TsilST* protos[16];
public:
  HH()
  {
  }


  HH(long double,long double,long double,long double,long double);
  
  HH(SMinput, long double);
  
  void init(long double,long double,long double,long double,long double);


  
  std::complex<long double> Mfin1(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParameters (&result, x, y, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMH); 
    return TSIL_GetFunction(&result, "M");
  }

  std::complex<long double> Vfin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMH); 
    return TSIL_GetBoldFunction(&result, "Vxzuv",0);
  }

  std::complex<long double> Ufin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMH); 
    return TSIL_GetBoldFunction(&result, "Uxzuv",0);
  }

  std::complex<long double> Tfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mu2); 
    TSIL_Evaluate( &result, MMH); 
    return TSIL_GetBoldFunction(&result, "Txuv",0);
  }

  std::complex<long double> Sfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mu2); 
    TSIL_Evaluate( &result, MMH); 
    return TSIL_GetBoldFunction(&result, "Sxuv",0);
  }

  
  void test()
  {
    std::vector<std::complex<long double> > diffMfin;
    std::vector<std::complex<long double> > diffVfin;
    std::vector<std::complex<long double> > diffUfin;
    std::vector<std::complex<long double> > diffTfin;
    std::vector<std::complex<long double> > diffSfin;
    std::vector<std::complex<long double> > diffIfin;


#include "testHH.hpp"
    
    for(int i = 0; i < diffMfin.size(); i++)
      std::cout << "Test diffMfin[" << i << "]= " << diffMfin[i] << std::endl;


    for(int i = 0; i < diffVfin.size(); i++)
      std::cout << "Test diffVfin[" << i << "]= " << diffVfin[i] << std::endl;


    for(int i = 0; i < diffUfin.size(); i++)
      std::cout << "Test diffUfin[" << i << "]= " << diffUfin[i] << std::endl;


    for(int i = 0; i < diffTfin.size(); i++)
      std::cout << "Test diffTfin[" << i << "]= " << diffTfin[i] << std::endl;


    for(int i = 0; i < diffSfin.size(); i++)
      std::cout << "Test diffSfin[" << i << "]= " << diffSfin[i] << std::endl;


  }
  std::complex<long double> m01(size_t nL = 2, size_t nH = 1);

  
  std::complex<long double> m10(size_t nL = 2, size_t nH = 1);

  
  std::complex<long double> m11(size_t nL = 2, size_t nH = 1);
  

  std::complex<long double> m20(size_t nL = 2, size_t nH = 1);

  // inline std::complex<long double> dd(size_t nG)
  // {
  //   return std::complex<long double>(nG == 0 ? 0 : 1,0);
  // }


  std::complex<long double> dalpha(long double MMt,long double MMH)
  {
    return 0;
  }
  
};


#endif  //  __HH_HPP__
