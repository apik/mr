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
bb::mygl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[14], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMW,-1);
    aryubbGL[4]=Tsil::I2(0,0,MMt,mu2);
    aryubbGL[5]=Tsil::A(MMt,mu2);
    aryubbGL[6]=pow(MMt,-1);
    aryubbGL[7]=Tsil::A(MMb,mu2);
    aryubbGL[8]=pow(MMb,-1);
    aryubbGL[9]=Tsil::Aeps(MMt,mu2);
    aryubbGL[10]=Tsil::Aeps(MMb,mu2);
   aryubbGL[11]=aryubbGL[5]*aryubbGL[6];
   aryubbGL[11]=9 - 7*aryubbGL[11];
   aryubbGL[11]=aryubbGL[5]*aryubbGL[11];
   aryubbGL[12]= - 3*aryubbGL[5] - 5./2.*MMt + 1./2.*MMH - 5*
   aryubbGL[7];
   aryubbGL[12]=aryubbGL[8]*aryubbGL[7]*aryubbGL[12];
   aryubbGL[11]=aryubbGL[11] + aryubbGL[12];
   aryubbGL[12]=aryubbGL[4] - aryubbGL[9];
   aryubbGL[13]=pow(Pi,2);
   aryubbGL[13]= - 5./2. - aryubbGL[13];
   aryubbGL[13]=MMt*aryubbGL[13];
   aryubbGL[11]=5*aryubbGL[10] + 1./3.*aryubbGL[13] - 1./12.*MMH + 1./2.
   *aryubbGL[11] - 4*aryubbGL[12];

      yubbGLret = aryubbGL[11]*aryubbGL[3]*pow(aryubbGL[2],2)*
      aryubbGL[1];
      return yubbGLret;
}
