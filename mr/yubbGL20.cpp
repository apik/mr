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
bb::ygl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubbGL[28], yubbGLret;

    aryubbGL[1]=double(boson);
    aryubbGL[2]=pow(SW,-1);
    aryubbGL[3]=pow(MMH,-1);
    aryubbGL[4]=pow(MMW,-1);
    aryubbGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryubbGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryubbGL[7]=pow(MMt,-1);
    aryubbGL[8]=Tsil::I2(0,MMH,MMt,mu2);
    aryubbGL[9]=Tsil::I2(0,0,MMH,mu2);
    aryubbGL[10]=Tsil::I2(0,0,MMt,mu2);
    aryubbGL[11]=Tsil::B(MMH,MMH,MMH,mu2);
    aryubbGL[12]=Tsil::B(MMH,MMt,MMt,mu2);
    aryubbGL[13]=Tsil::A(MMt,mu2);
    aryubbGL[14]=Tsil::B(MMt,MMt,MMH,mu2);
    aryubbGL[15]=Tsil::A(MMH,mu2);
    aryubbGL[16]=std::real(Tsil::B(0,0,MMt,mu2));
    aryubbGL[17]=Tsil::Aeps(MMH,mu2);
    aryubbGL[18]=Tsil::Aeps(MMt,mu2);
    aryubbGL[19]=std::real(Tsil::B(0,0,MMH,mu2));
   aryubbGL[20]=11*aryubbGL[12];
   aryubbGL[21]=pow(Pi,2);
   aryubbGL[22]= - 149./4. - 35*aryubbGL[21];
   aryubbGL[22]= - 11./4.*aryubbGL[16] + 1./8.*aryubbGL[22] - 
   aryubbGL[20];
   aryubbGL[23]=MMt*aryubbGL[3];
   aryubbGL[24]= - 3 - 1./2.*aryubbGL[21];
   aryubbGL[24]=aryubbGL[24]*aryubbGL[23];
   aryubbGL[25]=aryubbGL[3]*pow(aryubbGL[15],2);
   aryubbGL[25]= - 63./8.*aryubbGL[25] - 15./4.*aryubbGL[13] - 33*
   aryubbGL[18] - 13*aryubbGL[17];
   aryubbGL[25]=aryubbGL[3]*aryubbGL[25];
   aryubbGL[22]=29./8.*aryubbGL[24] + 1./4.*aryubbGL[22] + aryubbGL[25]
   ;
   aryubbGL[22]=MMt*aryubbGL[22];
   aryubbGL[24]=3*aryubbGL[12];
   aryubbGL[25]= - 3./4.*aryubbGL[16] + 85./16. - aryubbGL[24];
   aryubbGL[25]=aryubbGL[13]*aryubbGL[25];
   aryubbGL[26]= - 31./8. + 13*aryubbGL[23];
   aryubbGL[26]=aryubbGL[10]*aryubbGL[26];
   aryubbGL[25]=aryubbGL[25] + aryubbGL[26];
   aryubbGL[26]=pow(aryubbGL[13],2);
   aryubbGL[27]= - 11*aryubbGL[13] + 45./4.*aryubbGL[15];
   aryubbGL[27]=aryubbGL[15]*aryubbGL[27];
   aryubbGL[27]= - 51*aryubbGL[26] + aryubbGL[27];
   aryubbGL[27]=aryubbGL[3]*aryubbGL[27];
   aryubbGL[22]=aryubbGL[22] + 1./4.*aryubbGL[27] - 7./16.*aryubbGL[15]
    + 75./16.*aryubbGL[17] + 11*aryubbGL[18] - 17./8.*aryubbGL[6] + 1./
   2.*aryubbGL[25];
   aryubbGL[22]=MMt*aryubbGL[22];
   aryubbGL[20]=aryubbGL[20] + 197./4. + 9*aryubbGL[21];
   aryubbGL[20]=MMt*aryubbGL[20];
   aryubbGL[20]=aryubbGL[8] + aryubbGL[20] + 7*aryubbGL[18] + 5*
   aryubbGL[6];
   aryubbGL[25]=aryubbGL[13]*aryubbGL[7];
   aryubbGL[24]=7./2.*aryubbGL[25] - 13./4. + aryubbGL[24];
   aryubbGL[24]=aryubbGL[13]*aryubbGL[24];
   aryubbGL[25]=61./2. - aryubbGL[25];
   aryubbGL[25]=aryubbGL[15]*aryubbGL[25];
   aryubbGL[24]=aryubbGL[24] + aryubbGL[25];
   aryubbGL[25]=aryubbGL[18] + aryubbGL[6];
   aryubbGL[25]= - 3./2.*aryubbGL[8] + aryubbGL[9] + 1./2.*aryubbGL[25]
   ;
   aryubbGL[25]=aryubbGL[7]*aryubbGL[25];
   aryubbGL[21]=9./8.*aryubbGL[11] + 3./8.*aryubbGL[19] + 243./4.*S2 - 
   515./32. - 1./3.*aryubbGL[21] + aryubbGL[25];
   aryubbGL[25]=1./2.*MMH;
   aryubbGL[21]=aryubbGL[21]*aryubbGL[25];
   aryubbGL[20]=aryubbGL[21] - 15./8.*aryubbGL[5] - 17./8.*aryubbGL[9]
    + 7*aryubbGL[17] + 1./4.*aryubbGL[24] + 1./8.*aryubbGL[20];
   aryubbGL[20]=aryubbGL[20]*aryubbGL[25];
   aryubbGL[21]=7*aryubbGL[13] - 9*aryubbGL[15];
   aryubbGL[21]=aryubbGL[15]*aryubbGL[21];
   aryubbGL[21]=41./2.*aryubbGL[26] + aryubbGL[21];
   aryubbGL[24]=1./4.*MMt;
   aryubbGL[23]=49./4. - 11*aryubbGL[23];
   aryubbGL[23]=aryubbGL[8]*aryubbGL[23]*aryubbGL[24];
   aryubbGL[25]=pow(MMt,2);
   aryubbGL[24]=MMH*aryubbGL[24];
   aryubbGL[24]= - aryubbGL[25] + aryubbGL[24];
   aryubbGL[24]=aryubbGL[14]*aryubbGL[24];
   aryubbGL[20]=3./2.*aryubbGL[24] + aryubbGL[20] + aryubbGL[23] + 3./
   16.*aryubbGL[21] + aryubbGL[22];

      yubbGLret = 1./4.*aryubbGL[20]*pow(aryubbGL[4],2)*pow(
      aryubbGL[2],4)*aryubbGL[1];
      return yubbGLret;
}
