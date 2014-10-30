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
tt::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryutt[26], yuttret;

    aryutt[1]=double(nH);
    aryutt[2]=pow(CW,-1);
    aryutt[3]=pow(MMZ,-1);
    aryutt[4]=pow(SW,-1);
    aryutt[5]=Tsil::A(MMt,mu2);
    aryutt[6]=double(boson);
    aryutt[7]=Tsil::B(MMH,MMt,MMt,mu2);
    aryutt[8]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryutt[9]=pow(MMt,-1);
    aryutt[10]=Tsil::A(MMH,mu2);
    aryutt[11]=Tsil::A(MMZ,mu2);
    aryutt[12]=Tsil::A(MMW,mu2);
    aryutt[13]=std::real(Tsil::B(0,MMW,MMt,mu2));
    aryutt[14]=1/( - MMW + MMH);
   aryutt[15]=pow(aryutt[4],2);
   aryutt[16]=aryutt[15]*aryutt[12];
   aryutt[17]=7*aryutt[12];
   aryutt[18]=aryutt[17] + 3*aryutt[16];
   aryutt[18]=aryutt[18]*aryutt[15];
   aryutt[19]=pow(aryutt[2],2);
   aryutt[17]=aryutt[19]*aryutt[17];
   aryutt[20]=aryutt[19] + aryutt[15];
   aryutt[21]= - aryutt[10]*aryutt[20];
   aryutt[22]=aryutt[15] - 1;
   aryutt[23]= - aryutt[22]*aryutt[15];
   aryutt[23]=aryutt[23] + aryutt[19];
   aryutt[23]=aryutt[11]*aryutt[23];
   aryutt[17]=3*aryutt[23] + aryutt[21] + aryutt[18] + aryutt[17];
   aryutt[18]=1./4.*aryutt[13];
   aryutt[21]=aryutt[18] + aryutt[7];
   aryutt[21]=MMt*aryutt[20]*aryutt[21];
   aryutt[17]=1./4.*aryutt[17] + aryutt[21];
   aryutt[17]=aryutt[3]*aryutt[17];
   aryutt[21]=1./8.*aryutt[15];
   aryutt[23]=aryutt[21] + 17./72.*aryutt[19];
   aryutt[24]=aryutt[23] - 4./9.;
   aryutt[25]= - aryutt[8]*aryutt[24];
   aryutt[18]= - aryutt[22]*aryutt[18];
   aryutt[18]=aryutt[25] + aryutt[18];
   aryutt[18]=MMZ*aryutt[18];
   aryutt[22]= - aryutt[11]*aryutt[24];
   aryutt[16]=aryutt[18] - 1./4.*aryutt[16] + aryutt[22];
   aryutt[16]=aryutt[9]*aryutt[16];
   aryutt[18]= - aryutt[12] + aryutt[10];
   aryutt[18]=aryutt[14]*aryutt[18];
   aryutt[18]=3./8.*aryutt[18] + 1./8.*aryutt[13] - 3./16.;
   aryutt[15]=aryutt[15]*aryutt[18];
   aryutt[18]= - 7./72.*aryutt[19] + aryutt[21] + 8./9.;
   aryutt[18]=aryutt[8]*aryutt[18];
   aryutt[20]=aryutt[20]*aryutt[3];
   aryutt[21]=1./2. - aryutt[7];
   aryutt[21]=MMH*aryutt[21]*aryutt[20];
   aryutt[15]=1./8.*aryutt[21] + 1./2.*aryutt[17] + aryutt[18] + 7./144.
   *aryutt[19] - 8./9. + aryutt[15] + aryutt[16];
   aryutt[15]=aryutt[6]*aryutt[15];
   aryutt[16]=8./9. + aryutt[23];
   aryutt[16]=aryutt[9]*aryutt[16];
   aryutt[16]=aryutt[16] + 1./4.*aryutt[20];
   aryutt[16]=aryutt[6]*aryutt[16];
   aryutt[17]=aryutt[1]*aryutt[20];
   aryutt[16]= - 3./4.*aryutt[17] + aryutt[16];
   aryutt[16]=aryutt[5]*aryutt[16];
   aryutt[17]=MMt*aryutt[17];

      yuttret = aryutt[15] + aryutt[16] - 3./8.*aryutt[17];
      return yuttret;
}
