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

namespace mr
{
  template<class T>
  class WW 
  {
  };


  template<>
  class WW<OS> : public PoleMassAndCouplings
  {

    long double MMb, MMt, MMH, MMW, MMZ, mu2;
    long double SW, CW;

    std::unique_ptr<Tsil> protWHHWW;
    std::unique_ptr<Tsil> protWHZWW;
    std::unique_ptr<Tsil> protWZZWW;
    std::unique_ptr<Tsil> protWWHHH;
    std::unique_ptr<Tsil> protWWHZZ;
    std::unique_ptr<Tsil> protWWZZH;
    std::unique_ptr<Tsil> protWtZ00;
    std::unique_ptr<Tsil> protW0HWW;
    std::unique_ptr<Tsil> protW0Htt;
    std::unique_ptr<Tsil> protW0ZWW;
    std::unique_ptr<Tsil> protW0Ztt;
    std::unique_ptr<Tsil> protW0Z00;
    std::unique_ptr<Tsil> prot0WW0W;
    std::unique_ptr<Tsil> prot0Wt0t;
    std::unique_ptr<Tsil> prot0W0Z0;
    std::unique_ptr<Tsil> prot00Wt0;
    std::unique_ptr<Tsil> prot00W00;
    std::unique_ptr<Tsil> prot00ttZ;
    std::unique_ptr<Tsil> prot00tt0;
    std::unique_ptr<Tsil> prot0000Z;
    std::unique_ptr<Tsil> prot00000;
    std::unique_ptr<TsilSTU> protWH0H;
    std::unique_ptr<TsilSTU> protWZ0Z;
    std::unique_ptr<TsilSTU> protHW00;
  
  public:
    WW()
    {
    }

    WW(long double,long double,long double,long double,long double);

    WW(OSinput, long double);

    void init();
  

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
  class WW<MS> : public BaseMass
  {

    long double mmb, mmt, mmH, mmW, mmZ, mu2;
    long double s, c;

    std::unique_ptr<Tsil> protWHHWW;
    std::unique_ptr<Tsil> protWHZWW;
    std::unique_ptr<Tsil> protWZZWW;
    std::unique_ptr<Tsil> protWWHHH;
    std::unique_ptr<Tsil> protWWHZZ;
    std::unique_ptr<Tsil> protWWZZH;
    std::unique_ptr<Tsil> protWtZ00;
    std::unique_ptr<Tsil> protW0HWW;
    std::unique_ptr<Tsil> protW0Htt;
    std::unique_ptr<Tsil> protW0ZWW;
    std::unique_ptr<Tsil> protW0Ztt;
    std::unique_ptr<Tsil> protW0Z00;
    std::unique_ptr<Tsil> prot0WW0W;
    std::unique_ptr<Tsil> prot0Wt0t;
    std::unique_ptr<Tsil> prot0W0Z0;
    std::unique_ptr<Tsil> prot00Wt0;
    std::unique_ptr<Tsil> prot00W00;
    std::unique_ptr<Tsil> prot00ttZ;
    std::unique_ptr<Tsil> prot00tt0;
    std::unique_ptr<Tsil> prot0000Z;
    std::unique_ptr<Tsil> prot00000;
    std::unique_ptr<TsilSTU> protWH0H;
    std::unique_ptr<TsilSTU> protWZ0Z;
    std::unique_ptr<TsilSTU> protHW00;

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
} // namespace mr

#endif  //  __WW_HPP__
