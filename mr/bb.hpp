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

#ifndef __BB_HPP__
#define __BB_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"
#include "base.hpp"

class bb : public PoleMass
{

  long double MMb, MMt, MMH, MMW, MMZ, mu2;
  long double SW, CW;

  Tsil* prot0bb0b;

public:
  bb()
  {
  }

  bb(long double,long double,long double,long double,long double,long double);
  
  bb(OSinput, long double);
  
  void init(long double,long double,long double,long double,long double,long double);


  
  
  std::complex<long double> Mfin1(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;
    TSIL_SetParameters (&result, x, y, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMb); 
    return TSIL_GetFunction(&result, "M");
  }

  std::complex<long double> Vfin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMb); 
    return TSIL_GetBoldFunction(&result, "Vxzuv",0);
  }

  std::complex<long double> Ufin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMb); 
    return TSIL_GetBoldFunction(&result, "Uxzuv",0);
  }

  std::complex<long double> Tfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mu2); 
    TSIL_Evaluate( &result, MMb); 
    return TSIL_GetBoldFunction(&result, "Txuv",0);
  }

  std::complex<long double> Sfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mu2); 
    TSIL_Evaluate( &result, MMb); 
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


// #include "testbb.hpp"
    
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

  // 
  // Mass corrections
  // 
  long double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double x02(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double x03(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
  long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

  // Gaugeless limit
  long double xgl01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double xgl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  long double xgl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double xgl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


  // 
  // mass definition using Yukawa couplings 
  // mY=y/sqrt(2*sqrt(2)*GF)
  // 
  long double y01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double y02(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double y03(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double y10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double y11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  long double y20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


  // Gaugeless limit
  long double ygl01();

  long double ygl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double ygl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double ygl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


  std::complex<long double> det(const long double & a, const long double & b, const long double & c)
  {
    return 1./(a*a + b*b + c*c - 2*a*b - 2*b*c - 2*c*a);
  }
  
};


#endif  //  __BB_HPP__
