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

#include <bb.hpp>
std::complex<long double>
bb::ygl01(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[5], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=Tsil::A(mmb,mu2);
    aryubbGL[3]=pow(mmb,-1);
   aryubbGL[4]=aryubbGL[3]*aryubbGL[2];
   aryubbGL[4]= - 1./3. + aryubbGL[4];

      yubbGLret = 4*aryubbGL[4]*aryubbGL[1];
      return yubbGLret;
}
