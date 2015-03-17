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

#ifndef __TOP_HPP__
#define __TOP_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"
#include "base.hpp"

template<class T>
class tt
{
};

template<>
class tt<OS> : public PoleMassAndCouplings
{

long double MMb,MMt, MMH, MMW, MMZ, mu2;
long double SW, CW;

// alpha*alphaS
Tsil* protWt000;
Tsil* prot0ttHt;
Tsil* prot0ttZt;
Tsil* prot0tt0t;
TsilSTU* prottH0H;
TsilSTU* prottZ0Z;
  // alpha^2
  Tsil* protHHttH;
  Tsil* protHZttZ;
  Tsil* protHWt0W;
  Tsil* protHttHt;
  Tsil* protHttZt;
  Tsil* protZZttH;
  Tsil* protZWt0W;
  Tsil* protZttZt;
  Tsil* protZ0tW0;
  Tsil* protWW00Z;
  Tsil* protW00tW;
  Tsil* prot00WW0;
  TsilSTU* prot0W00;
  TsilST* prot000;

  // Gaugeless limit
  Tsil* protH0tt0;
  Tsil* protH0t00;
  Tsil* prot0Htt0;
  Tsil* prot0H0t0;
  Tsil* prot00ttH;
  Tsil* protHtt0t;
  Tsil* prot00t00;
  Tsil* prot000t0;

TsilST* protos[28];

// static const long double EPAIR2 = -1.; 

public:
tt()
{
}

  tt(long double,long double,long double,long double,long double);
  
  tt(OSinput, long double);
  
  void init();


  
  
  std::complex<long double> Mfin1(TSIL_REAL x, TSIL_REAL y, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;
    TSIL_SetParameters (&result, x, y, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetFunction(&result, "M");
  }

  std::complex<long double> Vfin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetBoldFunction(&result, "Vxzuv",0);
  }

  std::complex<long double> Ufin1(TSIL_REAL x, TSIL_REAL z, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersSTU (&result, x, z, u, v, mu2); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetBoldFunction(&result, "Uxzuv",0);
  }

  std::complex<long double> Tfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mu2); 
    TSIL_Evaluate( &result, MMt); 
    return TSIL_GetBoldFunction(&result, "Txuv",0);
  }

  std::complex<long double> Sfin1(TSIL_REAL x, TSIL_REAL u, TSIL_REAL v)
  {
    TSIL_DATA    result;

    TSIL_SetParametersST (&result, x, u, v, mu2); 
    TSIL_Evaluate( &result, MMt); 
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

#include "testtt.hpp"
// #include "testyt.hpp"
// #include "testttgl.hpp"
// #include "testytgl.hpp"

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

std::pair<long double,long double> test2(long double epsabs = 1.E-5,long double epsrel=0.001);
  // 
  // Pure QCD part m_ij=mY_ij by definition
  // 
  long double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double x02(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double x03(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  // 
  // Mass corrections
  // 
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
  long double ygl01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double ygl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double ygl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double ygl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
};


template<>
class tt<MS> : public PoleMass
{

long double mmb,mmt, mmH, mmW, mmZ, mu2;
long double s, c;

  // alpha*alphaS
  Tsil* protWt000;
  Tsil* prot0ttHt;
  Tsil* prot0ttZt;
  Tsil* prot0tt0t;
  TsilSTU* prottH0H;
  TsilSTU* prottZ0Z;
  // alpha^2
  Tsil* protHHttH;
  Tsil* protHZttZ;
  Tsil* protHWt0W;
  Tsil* protHttHt;
  Tsil* protHttZt;
  Tsil* protZZttH;
  Tsil* protZWt0W;
  Tsil* protZttZt;
  Tsil* protZ0tW0;
  Tsil* protWW00Z;
  Tsil* protW00tW;
  Tsil* prot00WW0;
  TsilSTU* prot0W00;
  TsilST* prot000;

  // Gaugeless limit
  Tsil* protH0tt0;
  Tsil* protH0t00;
  Tsil* prot0Htt0;
  Tsil* prot0H0t0;
  Tsil* prot00ttH;
  Tsil* protHtt0t;
  Tsil* prot00t00;
  Tsil* prot000t0;
  
  TsilST* protos[28];

// static const long double EPAIR2 = -1.; 

public:
  tt()
  {
  }
  
  tt(long double,long double,long double,long double,long double);
  
  tt(MSinput, long double);
  
  void init();


  
  
  // 
  // Pure QCD part m_ij=mY_ij by definition
  // 
  long double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double x02(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double x03(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  // 
  // Mass corrections
  // 
  long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

};


#endif  //  __TOP_HPP__
