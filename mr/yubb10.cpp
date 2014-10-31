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
bb::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubb[19], yubbret;

    aryubb[1]=double(boson);
    aryubb[2]=pow(CW,-1);
    aryubb[3]=pow(MMZ,-1);
    aryubb[4]=pow(SW,-1);
    aryubb[5]=Tsil::A(MMZ,mu2);
    aryubb[6]=Tsil::A(MMW,mu2);
    aryubb[7]=Tsil::A(MMt,mu2);
    aryubb[8]=Tsil::A(MMb,mu2);
    aryubb[9]=pow(MMb,-1);
    aryubb[10]=1/(MMt - MMW);
    aryubb[11]=1/( - MMW + MMH);
    aryubb[12]=Tsil::A(MMH,mu2);
   aryubb[13]=pow(aryubb[2],2);
   aryubb[14]= - 3*aryubb[6] - 1./4.*MMH + 5./4.*MMt + 3./2.*aryubb[7];
   aryubb[15]=5./6.*aryubb[5] - aryubb[14];
   aryubb[15]=aryubb[15]*aryubb[13];
   aryubb[16]=pow(aryubb[4],2);
   aryubb[17]=aryubb[6] - aryubb[5];
   aryubb[17]=aryubb[17]*aryubb[16];
   aryubb[17]=aryubb[17] + aryubb[5];
   aryubb[14]= - aryubb[14] + 3./2.*aryubb[17];
   aryubb[14]=aryubb[14]*aryubb[16];
   aryubb[14]=aryubb[15] + aryubb[14];
   aryubb[14]= - 1./3.*aryubb[5] + 1./4.*aryubb[14];
   aryubb[14]=aryubb[3]*aryubb[14];
   aryubb[15]=aryubb[12] - aryubb[6];
   aryubb[15]=aryubb[11]*aryubb[15];
   aryubb[17]=aryubb[6] - aryubb[7];
   aryubb[17]=aryubb[10]*aryubb[17];
   aryubb[17]=aryubb[17] + 1;
   aryubb[17]=MMZ*aryubb[17];
   aryubb[18]=aryubb[6] + aryubb[17];
   aryubb[18]=aryubb[10]*aryubb[18];
   aryubb[15]=aryubb[18] + aryubb[15];
   aryubb[15]=aryubb[15]*aryubb[16];
   aryubb[16]=aryubb[10]*aryubb[17];
   aryubb[15]=aryubb[16] - aryubb[15];
   aryubb[16]=aryubb[8]*aryubb[9];
   aryubb[16]= - 1./2. + aryubb[16];
   aryubb[13]=aryubb[14] - 11./72.*aryubb[13] + 1./3.*aryubb[16] - 3./8.
   *aryubb[15];

      yubbret = aryubb[13]*aryubb[1];
      return yubbret;
}
