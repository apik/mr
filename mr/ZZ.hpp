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

#ifndef __ZZ_HPP__
#define __ZZ_HPP__

#include "tsil.hpp"
#include "sminput.hpp"
#include "operators.hpp"
#include "constants.hpp"
#include "base.hpp"

namespace mr
{
  template<class T>
  class ZZ 
  {
  };

  template<>
  class ZZ<OS> : public PoleMassAndCouplings
  {

    long double MMb, MMt, MMH, MMW, MMZ, mu2;
    long double SW, CW;

    std::unique_ptr<Tsil> protZHHZZ;
    std::unique_ptr<Tsil> protZZHHH;
    std::unique_ptr<Tsil> protZWHWW;
    std::unique_ptr<Tsil> prottZtHt;
    std::unique_ptr<Tsil> protWWWWH;
    std::unique_ptr<Tsil> protWWWWZ;
    std::unique_ptr<Tsil> protWWWW0;
    std::unique_ptr<Tsil> protWtWt0;
    std::unique_ptr<Tsil> protW0W0t;
    std::unique_ptr<Tsil> protW0W00;
    std::unique_ptr<Tsil> protttttH;
    std::unique_ptr<Tsil> protttttZ;
    std::unique_ptr<Tsil> prottttt0;
    std::unique_ptr<Tsil> prott0t0W;
    std::unique_ptr<Tsil> prot0000Z;
    std::unique_ptr<Tsil> prot0000W;
    std::unique_ptr<Tsil> prot00000;
    std::unique_ptr<Tsil> protWZWHW;
    std::unique_ptr<TsilSTU> protHZ00;


  public:
    ZZ()
    {
    }

    ZZ(long double,long double,long double,long double,long double);

    ZZ(OSinput, long double);

    void init();
  
    long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
    long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

    long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


    long double y10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
    long double y11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

    long double y20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
    
  };

  template<>
  class ZZ<MS> : public BaseMass
  {

    long double mmb, mmt, mmH, mmW, mmZ, mu2;
    long double s, c;

    std::unique_ptr<Tsil> protZHHZZ;
    std::unique_ptr<Tsil> protZZHHH;
    std::unique_ptr<Tsil> protZWHWW;
    std::unique_ptr<Tsil> prottZtHt;
    std::unique_ptr<Tsil> protWWWWH;
    std::unique_ptr<Tsil> protWWWWZ;
    std::unique_ptr<Tsil> protWWWW0;
    std::unique_ptr<Tsil> protWtWt0;
    std::unique_ptr<Tsil> protW0W0t;
    std::unique_ptr<Tsil> protW0W00;
    std::unique_ptr<Tsil> protttttH;
    std::unique_ptr<Tsil> protttttZ;
    std::unique_ptr<Tsil> prottttt0;
    std::unique_ptr<Tsil> prott0t0W;
    std::unique_ptr<Tsil> prot0000Z;
    std::unique_ptr<Tsil> prot0000W;
    std::unique_ptr<Tsil> prot00000;
    std::unique_ptr<Tsil> protWZWHW;
    std::unique_ptr<TsilSTU> protHZ00;
  
  public:
    ZZ()
    {
    }

    ZZ(long double,long double,long double,long double,long double);

    ZZ(MSinput, long double);

    void init();


    long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
    long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
    

    long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    
  };
} // namespace mr

#endif  //  __ZZ_HPP__
