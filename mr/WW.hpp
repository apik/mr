//
// MR - 2-loop matching and 3-loop Running, including full 2-loop EW corrections
// Copyright (C) 2014 Andrey Pikelner <pikelner@theor.jinr.ru>
//
// This file is part of MR.
//
// MR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MR.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __WW_HPP__
#define __WW_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"
#include "base.hpp"


template<class T>
class WW 
{
};


template<>
class WW<OS> : public PoleMassAndCouplings
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
  TsilSTU* protWH0H;
  TsilSTU* protWZ0Z;
  TsilSTU* protHW00;
  
  TsilST* protos[24];
public:
  WW()
  {
  }

  WW(long double,long double,long double,long double,long double);

  WW(OSinput, long double);

  void init();
  
  std::complex<long double> Mfin1(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParameters (&result, x, y, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMW); 
    return 0*TSIL_GetFunction(&result, "M");
  }

  std::complex<long double> Vfin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMW); 
    return 0*TSIL_GetBoldFunction(&result, "Vxzuv",0);
  }

  std::complex<long double> Ufin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMW); 
    return 0*TSIL_GetBoldFunction(&result, "Uxzuv",0);
  }

  std::complex<long double> Tfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mu2); 
    TSIL_Evaluate( &result, MMW); 
    return 0*TSIL_GetBoldFunction(&result, "Txuv",0);
  }

  std::complex<long double> Sfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mu2); 
    TSIL_Evaluate( &result, MMW); 
    return 0*TSIL_GetBoldFunction(&result, "Sxuv",0);
  }


  std::pair<long double,long double> test2(long double epsabs = 1.E-5,long double epsrel=0.001);  

  void test()
  {
    std::vector<std::complex<long double> > diffMfin;
    std::vector<std::complex<long double> > diffVfin;
    std::vector<std::complex<long double> > diffUfin;
    std::vector<std::complex<long double> > diffTfin;
    std::vector<std::complex<long double> > diffSfin;
    std::vector<std::complex<long double> > diffIfin;


// #include "testWW.hpp"
    
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

  // Pole -> MS
  
  long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

  long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


  // Yukawa
  long double y10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  long double y11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

  long double y20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

};


template<>
class WW<MS>  : public PoleMass 
{

  long double mmb, mmt, mmH, mmW, mmZ, mu2;
  long double s, c;

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
  TsilSTU* protWH0H;
  TsilSTU* protWZ0Z;
  TsilSTU* protHW00;
  
  TsilST* protos[24];
public:
  WW()
  {
  }

  WW(long double,long double,long double,long double,long double);

  WW(MSinput, long double);

  void init();
  
  // Pole -> MS
  
  long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

  long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


};


#endif  //  __WW_HPP__
