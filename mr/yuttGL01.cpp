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
tt::mygl01(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[5], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=Tsil::A(MMt,mu2);
    aryuttGL[3]=pow(MMt,-1);
   aryuttGL[4]=aryuttGL[2]*aryuttGL[3];
   aryuttGL[4]= - 1./3. + aryuttGL[4];

      yuttGLret = 4*aryuttGL[4]*aryuttGL[1];
      return yuttGLret;
}
