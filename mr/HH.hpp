////
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

namespace mr
{
  template<class T>
  class HH 
  {
  };

  template<>
  class HH<OS> : public PoleMassAndCouplings
  {

    double MMb, MMt, MMH, MMW, MMZ, mu2;
    double SW, CW;
  
    std::unique_ptr<Tsil> protHHHHH;
    std::unique_ptr<Tsil> protHZHZZ;
    std::unique_ptr<Tsil> protHWHWW;
    std::unique_ptr<Tsil> protHtHtt;
    std::unique_ptr<Tsil> protZZZZH;
    std::unique_ptr<Tsil> protZWZWW;
    std::unique_ptr<Tsil> protZtZtt;
    std::unique_ptr<Tsil> protWWWWH;
    std::unique_ptr<Tsil> protWWWWZ;
    std::unique_ptr<Tsil> protWWWW0;
    std::unique_ptr<Tsil> protWtWt0;
    std::unique_ptr<Tsil> protttttH;
    std::unique_ptr<Tsil> protttttZ;
    std::unique_ptr<Tsil> prottttt0;
    std::unique_ptr<TsilSTU> protZZ00;
    std::unique_ptr<TsilSTU> protWW00;
  
    // Gaugeless limit
    std::unique_ptr<Tsil> prot0H0H0;
    std::unique_ptr<Tsil> prot0t0tt;
    std::unique_ptr<Tsil> prot0t0t0;
    std::unique_ptr<Tsil> prot0000H;

  public:
    HH()
    {
    }


    HH(double,double,double,double,double);
  
    HH(OSinput, double);
  
    void init();

    double Q()
    {
      return sqrt(mu2);
    }

    double QQ()
    {
      return mu2;
    }
  
  
    double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
    double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

    double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


    // Gaugeless limit
    double xgl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double xgl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double xgl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


    
    // 
    // \lambda/sqrt(2)/GF/MMH
    // 
    double y10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double y11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double y20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    // Gaugeless limit

    double ygl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double ygl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double ygl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
  };


  template<>
  class HH<MS> : public BaseMass
  {

    double mmb, mmt, mmH, mmW, mmZ, mu2;
    double s, c;
  
    std::unique_ptr<Tsil> protHHHHH;
    std::unique_ptr<Tsil> protHZHZZ;
    std::unique_ptr<Tsil> protHWHWW;
    std::unique_ptr<Tsil> protHtHtt;
    std::unique_ptr<Tsil> protZZZZH;
    std::unique_ptr<Tsil> protZWZWW;
    std::unique_ptr<Tsil> protZtZtt;
    std::unique_ptr<Tsil> protWWWWH;
    std::unique_ptr<Tsil> protWWWWZ;
    std::unique_ptr<Tsil> protWWWW0;
    std::unique_ptr<Tsil> protWtWt0;
    std::unique_ptr<Tsil> protttttH;
    std::unique_ptr<Tsil> protttttZ;
    std::unique_ptr<Tsil> prottttt0;
    std::unique_ptr<TsilSTU> protZZ00;
    std::unique_ptr<TsilSTU> protWW00;
  
    // Gaugeless limit
    std::unique_ptr<Tsil> prot0H0H0;
    std::unique_ptr<Tsil> prot0t0tt;
    std::unique_ptr<Tsil> prot0t0t0;
    std::unique_ptr<Tsil> prot0000H;

  public:
    HH()
    {
    }


    HH(double,double,double,double,double);
  
    HH(MSinput, double);
  
    void init();

    //  double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
    double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
    double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

    double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


  };

} // namespace mr

#endif  //  __HH_HPP__
