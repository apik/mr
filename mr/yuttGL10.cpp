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
tt::ygl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[10], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=pow(SW,-1);
    aryuttGL[3]=pow(MMW,-1);
    aryuttGL[4]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[5]=Tsil::A(MMH,mu2);
    aryuttGL[6]=Tsil::A(MMt,mu2);
    aryuttGL[7]=std::real(Tsil::B(0,0,MMt,mu2));
   aryuttGL[8]=1./2. - aryuttGL[4];
   aryuttGL[8]=MMH*aryuttGL[8];
   aryuttGL[8]= - aryuttGL[5] + aryuttGL[8];
   aryuttGL[9]= - 3 + aryuttGL[7];
   aryuttGL[9]=1./4.*aryuttGL[9] + aryuttGL[4];
   aryuttGL[9]=MMt*aryuttGL[9];
   aryuttGL[8]=aryuttGL[9] + 1./4.*aryuttGL[8] - aryuttGL[6];

      yuttGLret = 1./2.*aryuttGL[8]*aryuttGL[3]*pow(aryuttGL[2],2)*
      aryuttGL[1];
      return yuttGLret;
}
