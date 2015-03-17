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

#ifndef __HH_HPP__
#define __HH_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"
#include "base.hpp"

template<class T>
class HH 
//: public PoleMass
{
};

template<>
class HH<OS> : public PoleMassAndCouplings
{

  long double MMb, MMt, MMH, MMW, MMZ, mu2;
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
  
  // Gaugeless limit
  Tsil* prot0H0H0;
  Tsil* prot0t0tt;
  Tsil* prot0t0t0;
  Tsil* prot0000H;

  TsilST* protos[20];
public:
  HH()
  {
  }


  HH(long double,long double,long double,long double,long double);
  
  HH(OSinput, long double);
  
  void init();

  long double Q()
  {
    return sqrt(mu2);
  }

  long double QQ()
  {
    return mu2;
  }
  
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
  
  long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

  long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


  // Gaugeless limit
  long double xgl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  long double xgl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double xgl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


    
  // 
  // \lambda/sqrt(2)/GF/MMH
  // 
  long double y10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double y11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  long double y20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  // Gaugeless limit

  long double ygl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double ygl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double ygl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
};


template<>
class HH<MS> : public PoleMass
{

  long double mmb, mmt, mmH, mmW, mmZ, mu2;
  long double s, c;
  
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
  
  // Gaugeless limit
  Tsil* prot0H0H0;
  Tsil* prot0t0tt;
  Tsil* prot0t0t0;
  Tsil* prot0000H;

  TsilST* protos[20];
public:
  HH()
  {
  }


  HH(long double,long double,long double,long double,long double);
  
  HH(MSinput, long double);
  
  void init();

//  long double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

  long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


  // Gaugeless limit
  long double xgl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  long double xgl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double xgl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
};



#endif  //  __HH_HPP__
