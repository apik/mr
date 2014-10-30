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

#include <tt.hpp>
std::complex<long double>
tt::mgl01(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[5], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=Tsil::A(MMt,mu2);
    armttbarGL[3]=pow(MMt,-1);
   armttbarGL[4]=armttbarGL[2]*armttbarGL[3];
   armttbarGL[4]= - 1./3. + armttbarGL[4];

      mttbarGLret = 4*armttbarGL[4]*armttbarGL[1];
      return mttbarGLret;
}
