#ifndef __TSIL_HPP__
#define __TSIL_HPP__
// #include <stdio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>



#ifdef __cplusplus
extern "C"{
#endif  

#include "tsil.h"  

#ifdef __cplusplus
}
#endif

// workaround for using C99 complex with C++ std::complex
#ifdef complex                 
#undef complex
#endif

// Workaround using gcc extension to work both in gcc and clang
// http://clang.llvm.org/docs/LanguageExtensions.html#initializer-lists-for-complex-numbers-in-c
std::complex<long double> c2pp(TSIL_COMPLEX);

class TsilST
{
protected:
  TSIL_DATA    result;
  bool      evaluated;
  TSIL_REAL     scale;
  TSIL_REAL         s;
public:
  
  TsilST()
  {
    evaluated = false;
  }

  TsilST(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq)
  {
    TSIL_SetParametersST (&result, x, y, z, qq); 
    evaluated = false;
  }

  TsilST(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL qq, TSIL_REAL s_): scale(qq), s(s_)
  {
    TSIL_SetParametersST (&result, x, y, z, qq); 
    this->evaluate(s); 
  }

  void evaluate(TSIL_REAL s)
  {
    TSIL_Evaluate (&result, s);
    evaluated = true;
  }

  static std::complex<long double> A(TSIL_REAL m, TSIL_REAL qq)
  {
    return c2pp(TSIL_A(m,qq));
  }

  static std::complex<long double> Aeps(TSIL_REAL m, TSIL_REAL qq)
  {
    return c2pp(TSIL_Aeps(m,qq));
  }

  static std::complex<long double> B(TSIL_REAL m1, TSIL_REAL m2,TSIL_REAL s,TSIL_REAL qq)
  {
    return c2pp(TSIL_B(m1, m2, s, qq));
  }

  static std::complex<long double> Beps(TSIL_REAL m1, TSIL_REAL m2,TSIL_REAL s,TSIL_REAL qq)
  {
    return c2pp(TSIL_Beps(m1, m2, s, qq));
  }

  static std::complex<long double> I2(TSIL_REAL x, TSIL_REAL y,TSIL_REAL z,TSIL_REAL qq)
  {
    return c2pp(TSIL_Aeps(x, qq)) + c2pp(TSIL_Aeps(y, qq)) + c2pp(TSIL_Aeps(z, qq)) + c2pp(TSIL_I2(x, y, z, qq));
  }

  // 
  // S-type
  // 
  inline std::complex<long double> Svyz(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Svyz", -epsord));
  }

  inline std::complex<long double> Suxv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Suxv", -epsord));
  }

  // 
  // T-type
  // 
  inline std::complex<long double> Tvyz(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Tvyz", -epsord));
  }

  inline std::complex<long double> Tuxv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Tuxv", -epsord));
  }

  inline std::complex<long double> Tyzv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Tyzv", -epsord));
  }
 
  inline std::complex<long double> Txuv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Txuv", -epsord));
  }
  
  inline std::complex<long double> Tzyv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Tzyv", -epsord));
  }

  inline std::complex<long double> Tvxu(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Tvxu", -epsord));
  }

  // GiNaC version!!!!!

  // static GiNaC::ex Aginac(TSIL_REAL m, TSIL_REAL qq)
  // {
  //   return g2pp(TSIL_A(m,qq));
  // }

  // static GiNaC::ex Aepsginac(TSIL_REAL m, TSIL_REAL qq)
  // {
  //   return g2pp(TSIL_Aeps(m,qq));
  // }

  // static GiNaC::ex Bginac(TSIL_REAL m1, TSIL_REAL m2,TSIL_REAL s,TSIL_REAL qq)
  // {
  //   return g2pp(TSIL_B(m1, m2, s, qq));
  // }

  // static GiNaC::ex Bepsginac(TSIL_REAL m1, TSIL_REAL m2,TSIL_REAL s,TSIL_REAL qq)
  // {
  //   return g2pp(TSIL_Beps(m1, m2, s, qq));
  // }

  // static GiNaC::ex I2ginac(TSIL_REAL x, TSIL_REAL y,TSIL_REAL z,TSIL_REAL qq)
  // {
  //   return g2pp(TSIL_Aeps(x, qq)) + g2pp(TSIL_Aeps(y, qq)) + g2pp(TSIL_Aeps(z, qq)) + g2pp(TSIL_I2(x, y, z, qq));
  // }

  // // 
  // // S-type
  // // 
  // inline std::complex<long double> Svyz(int epsord)
  // {
  //   return c2pp(TSIL_GetBoldFunction(&result, "Svyz", -epsord));
  // }

  // inline GiNaC::ex Suxvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Suxv", -epsord));
  // }

  // 
  // T-type
  // 
  // inline GiNaC::ex Tvyzginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Tvyz", -epsord));
  // }

  // inline GiNaC::ex Tuxvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Tuxv", -epsord));
  // }

  // inline GiNaC::ex Tyzvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Tyzv", -epsord));
  // }
 
  // inline GiNaC::ex Txuvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Txuv", -epsord));
  // }
  
  // inline GiNaC::ex Tzyvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Tzyv", -epsord));
  // }

  // inline GiNaC::ex Tvxuginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Tvxu", -epsord));
  // }


  void print()
  {
    if(evaluated) TSIL_PrintData (&result);
    else std::cout << "Integral unevaluated!!!" <<  std::endl;
  }

};


class TsilSTU : public TsilST 
{
  public:

  TsilSTU()
  {
    evaluated = false;
  }

  TsilSTU(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v, TSIL_REAL qq)
  {
    TSIL_SetParametersSTU (&result, x, z, u, v, qq); 
    evaluated = false;
  }

  TsilSTU(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v, TSIL_REAL qq, TSIL_REAL s_)
  {
    scale = qq;
    s = s_;
    TSIL_SetParametersSTU (&result, x, z, u, v, qq); 
    this->evaluate(s); 
  }

  // 
  // U-type
  // 
  inline std::complex<long double> Uzxyv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Uzxyv", -epsord));
  }

  inline std::complex<long double> Uuyxv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Uuyxv", -epsord));
  }

  inline std::complex<long double> Uxzuv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Uxzuv", -epsord));
  }

  inline std::complex<long double> Uyuzv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Uyuzv", -epsord));
  }


  // Ginac

  // inline GiNaC::ex Uzxyvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Uzxyv", -epsord));
  // }

  // inline GiNaC::ex Uuyxvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Uuyxv", -epsord));
  // }

  // inline GiNaC::ex Uxzuvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Uxzuv", -epsord));
  // }

  // inline GiNaC::ex Uyuzvginac(int epsord)
  // {
  //   return g2pp(TSIL_GetBoldFunction(&result, "Uyuzv", -epsord));
  // }

  
  // 
  // V-type
  // 
  inline std::complex<long double> Vzxyv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Vzxyv", -epsord));
  }

  inline std::complex<long double> Vuyxv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Vuyxv", -epsord));
  }

  inline std::complex<long double> Vxzuv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Vxzuv", -epsord));
  }

  inline std::complex<long double> Vyuzv(int epsord)
  {
    return c2pp(TSIL_GetBoldFunction(&result, "Vyuzv", -epsord));
  }


};

class Tsil : public TsilSTU
{

public:
  Tsil(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v, TSIL_REAL qq)
  {
    TSIL_SetParameters (&result, x, y, z, u, v, qq); 
    evaluated = false;
  }

  Tsil(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v, TSIL_REAL qq, TSIL_REAL s_)
  {
    scale = qq;
    s = s_;  
    TSIL_SetParameters (&result, x, y, z, u, v, qq); 
    this->evaluate(s); 
  }

  void evaluate(TSIL_REAL s)
  {
    TSIL_Evaluate (&result, s);
    evaluated = true;
  }

  // 
  // M-type
  // 
  inline std::complex<long double> M(int epsord=0)
  {
    return c2pp(TSIL_GetFunction(&result, "M"));
  }

  // // GiNaC
  // inline GiNaC::ex Mginac(int epsord=0)
  // {
  //   return g2pp(TSIL_GetFunction(&result, "M"));
  // }



  void print()
  {
    if(evaluated) TSIL_PrintData (&result);
    else std::cout << "Integral unevaluated!!!" <<  std::endl;
  }

};



std::complex<long double> csqrt(long double);
std::complex<long double> Li2(std::complex<long double>);
std::complex<long double> Li3(std::complex<long double>);

std::complex<long double> acc(std::complex<long double>);
std::complex<long double> inv(std::complex<long double>);




#endif  // __TSIL_HPP__
