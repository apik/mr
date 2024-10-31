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

namespace mr
{
  template<class T>
  class bb
  {
  };

  template<>
  class bb<OS> : public PoleMassAndCouplings
  {

    double MMb, MMt, MMH, MMW, MMZ, mu2;
    double SW, CW;

    std::unique_ptr<Tsil> prot0bb0b;

  public:
    bb()
    {
    }

    bb(double,double,double,double,double,double);
  
    bb(OSinput, double);
  
    void init(double,double,double,double,double,double);

    // 
    // Mass corrections
    // 
    double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double x02(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double x03(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double x04(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  
    double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  

    // Gaugeless limit
    double xgl01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double xgl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double xgl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double xgl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


    // 
    // mass definition using Yukawa couplings 
    // mY=y/sqrt(2*sqrt(2)*GF)
    // 
    double y01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double y02(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double y03(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double y04(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double y10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double y11(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double y20(size_t nL = 2, size_t nH = 1, size_t boson = 1);


    // Gaugeless limit
    double ygl01();

    double ygl10(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double ygl11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double ygl20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  private:
    std::complex<double> det(const double & a, const double & b, const double & c)
    {
      return 1./(a*a + b*b + c*c - 2*a*b - 2*b*c - 2*c*a);
    }

  
  };


  template<>
  class bb<MS> : public BaseMass
  {

    double mmb,mmt, mmH, mmW, mmZ, mu2;
    double s, c;
  
    std::unique_ptr<Tsil> prot0bb0b;

  public:
    bb()
    {
    }

    bb(double,double,double,double,double,double);
  
    bb(MSinput, double);
  
    void init();

    // 
    // Pure QCD part m_ij=mY_ij by definition
    // 
    double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double x02(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double x03(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double x04(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    // 
    // Mass corrections
    // 
    double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1);
  
    double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1);

    double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1);

  private:  
    std::complex<double> det(const double & a, const double & b, const double & c)
    {
      return 1./(a*a + b*b + c*c - 2*a*b - 2*b*c - 2*c*a);
    }
  };
} // namespace mr

#endif  //  __BB_HPP__
