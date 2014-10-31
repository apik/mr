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

#ifndef __BASE_HPP__
#define __BASE_HPP__

class PoleMass
{
    // 
  // Pure QCD part m_ij=mY_ij by definition
  // 
  // virtual std::complex<long double> m01();

  // virtual std::complex<long double> m02(size_t nL = 5);

  // virtual std::complex<long double> m03(size_t nL = 5);

  // 
  // Mass corrections
  // 
public:
  virtual std::complex<long double> m10(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;
  
  virtual std::complex<long double> m11(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

  virtual std::complex<long double> m20(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

  virtual std::complex<long double> m02(size_t nL = 5)
  {
    std::cout << "Order a^0*as^2 is not implemented for this particle" << std::endl;
  }

  // and meta method
  long double m(size_t apow, size_t aspow, size_t nL = 2, size_t nH = 1, size_t boson = 1)
  {
    if(apow == 1 && aspow == 0)
      return m10(nL, nH, boson).real();
    if(apow == 1 && aspow == 1)
      return m11(nL, nH, boson).real();
    if(apow == 2 && aspow == 0)
      return m20(nL, nH, boson).real();

  }
  // Gaugeless limit
  // virtual std::complex<long double> mgl01(size_t nL = 2, size_t nH = 1);

  // virtual std::complex<long double> mgl10(size_t nL = 2, size_t nH = 1);
  
  // virtual std::complex<long double> mgl11(size_t nL = 2, size_t nH = 1);

  // virtual std::complex<long double> mgl20(size_t nL = 2, size_t nH = 1);
  

  // 
  // mass definition using Yukawa couplings 
  // mY=y/sqrt(2*sqrt(2)*GF)
  // 
  // virtual std::complex<long double> my01() = 0;

  virtual std::complex<long double> my10(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

  virtual std::complex<long double> my11(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

  virtual std::complex<long double> my20(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

  long double my(size_t apow, size_t aspow, size_t nL = 2, size_t nH = 1, size_t boson = 1)
  {
    if(apow == 1 && aspow == 0)
      return my10(nL, nH, boson).real();
    if(apow == 1 && aspow == 1)
      return my11(nL, nH, boson).real();
    if(apow == 2 && aspow == 0)
      return my20(nL, nH, boson).real();

  }

  // Gaugeless limit
  // virtual std::complex<long double> mygl01();

  // virtual std::complex<long double> mygl10(size_t nL = 2, size_t nH = 1);

  // virtual std::complex<long double> mygl11(size_t nL = 2, size_t nH = 1);

  // virtual std::complex<long double> mygl20(size_t nL = 2, size_t nH = 1);


};

// class RunningMass
// {
//   virtual void  m01()
// };

#endif  // __BASE_HPP__
