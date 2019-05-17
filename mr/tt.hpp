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


namespace mr
{
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
    std::unique_ptr<Tsil> protWt000;
    std::unique_ptr<Tsil> prot0ttHt;
    std::unique_ptr<Tsil> prot0ttZt;
    std::unique_ptr<Tsil> prot0tt0t;
    std::unique_ptr<TsilSTU> prottH0H;
    std::unique_ptr<TsilSTU> prottZ0Z;
    // alpha^2
    std::unique_ptr<Tsil> protHHttH;
    std::unique_ptr<Tsil> protHZttZ;
    std::unique_ptr<Tsil> protHWt0W;
    std::unique_ptr<Tsil> protHttHt;
    std::unique_ptr<Tsil> protHttZt;
    std::unique_ptr<Tsil> protZZttH;
    std::unique_ptr<Tsil> protZWt0W;
    std::unique_ptr<Tsil> protZttZt;
    std::unique_ptr<Tsil> protZ0tW0;
    std::unique_ptr<Tsil> protWW00Z;
    std::unique_ptr<Tsil> protW00tW;
    std::unique_ptr<Tsil> prot00WW0;
    std::unique_ptr<TsilSTU> prot0W00;
    std::unique_ptr<TsilST> prot000;

    // Gaugeless limit
    std::unique_ptr<Tsil> protH0tt0;
    std::unique_ptr<Tsil> protH0t00;
    std::unique_ptr<Tsil> prot0Htt0;
    std::unique_ptr<Tsil> prot0H0t0;
    std::unique_ptr<Tsil> prot00ttH;
    std::unique_ptr<Tsil> protHtt0t;
    std::unique_ptr<Tsil> prot00t00;
    std::unique_ptr<Tsil> prot000t0;

    // std::unique_ptr<TsilST> protos[28];
  
  public:
    tt()
    {
    }

    tt(long double,long double,long double,long double,long double);
  
    tt(OSinput, long double);
  
    void init();

    // 
    // Pure QCD part m_ij=mY_ij by definition
    // 
    long double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    long double x02(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    long double x03(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    long double x04(size_t nL = 2, size_t nH = 1, size_t boson = 1);

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

    long double y04(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    // EW
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
  class tt<MS> : public BaseMass
  {

    long double mmb,mmt, mmH, mmW, mmZ, mu2;
    long double s, c;

    // alpha*alphaS
    std::unique_ptr<Tsil> protWt000;
    std::unique_ptr<Tsil> prot0ttHt;
    std::unique_ptr<Tsil> prot0ttZt;
    std::unique_ptr<Tsil> prot0tt0t;
    std::unique_ptr<TsilSTU> prottH0H;
    std::unique_ptr<TsilSTU> prottZ0Z;
    // alpha^2
    std::unique_ptr<Tsil> protHHttH;
    std::unique_ptr<Tsil> protHZttZ;
    std::unique_ptr<Tsil> protHWt0W;
    std::unique_ptr<Tsil> protHttHt;
    std::unique_ptr<Tsil> protHttZt;
    std::unique_ptr<Tsil> protZZttH;
    std::unique_ptr<Tsil> protZWt0W;
    std::unique_ptr<Tsil> protZttZt;
    std::unique_ptr<Tsil> protZ0tW0;
    std::unique_ptr<Tsil> protWW00Z;
    std::unique_ptr<Tsil> protW00tW;
    std::unique_ptr<Tsil> prot00WW0;
    std::unique_ptr<TsilSTU> prot0W00;
    std::unique_ptr<TsilST> prot000;

    // Gaugeless limit
    std::unique_ptr<Tsil> protH0tt0;
    std::unique_ptr<Tsil> protH0t00;
    std::unique_ptr<Tsil> prot0Htt0;
    std::unique_ptr<Tsil> prot0H0t0;
    std::unique_ptr<Tsil> prot00ttH;
    std::unique_ptr<Tsil> protHtt0t;
    std::unique_ptr<Tsil> prot00t00;
    std::unique_ptr<Tsil> prot000t0;
  
    // std::unique_ptr<TsilST> protos[28];

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

    long double x04(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    // 
    // Mass corrections
    // 
    long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  };
} // namespace mr

#endif  //  __TOP_HPP__
