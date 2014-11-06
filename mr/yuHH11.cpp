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

#include <HH.hpp>
std::complex<long double>
HH<OS>::y11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHH[20], yuHHret;

    aryuHH[1]=double(nH);
    aryuHH[2]=pow(CW,-1);
    aryuHH[3]=pow(MMH,-1);
    aryuHH[4]=pow(MMZ,-1);
    aryuHH[5]=pow(SW,-1);
    aryuHH[6]=Tsil::I2(0,0,MMt,mu2);
    aryuHH[7]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHH[8]=Tsil::A(MMt,mu2);
    aryuHH[9]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHH[10]=pow(MMt,-1);
    aryuHH[11]=Tsil::Aeps(MMt,mu2);
    aryuHH[12]=prottttt0->M(0);
    aryuHH[13]=prottttt0->Vxzuv(0);
    aryuHH[14]=prottttt0->Suxv(0);
   aryuHH[15]=2*aryuHH[13];
   aryuHH[16]=aryuHH[15] + aryuHH[12];
   aryuHH[17]=MMt*aryuHH[3];
   aryuHH[18]=8*aryuHH[17];
   aryuHH[16]=aryuHH[16]*aryuHH[18];
   aryuHH[18]=10 + 7*aryuHH[7];
   aryuHH[18]=aryuHH[18]*aryuHH[7];
   aryuHH[18]=aryuHH[18] + 25 - 4*aryuHH[9];
   aryuHH[18]=aryuHH[18]*aryuHH[3];
   aryuHH[15]=aryuHH[15] + 3*aryuHH[12];
   aryuHH[15]=aryuHH[18] - aryuHH[16] + 2*aryuHH[15];
   aryuHH[16]=4*MMt;
   aryuHH[15]=aryuHH[15]*aryuHH[16];
   aryuHH[16]=4 + 3*aryuHH[7];
   aryuHH[18]=2*aryuHH[7];
   aryuHH[16]=aryuHH[16]*aryuHH[18];
   aryuHH[18]=MMH*aryuHH[12];
   aryuHH[19]=aryuHH[3]*aryuHH[11];
   aryuHH[15]=aryuHH[15] - 4*aryuHH[18] - aryuHH[16] - 85./2. - 24*
   aryuHH[19];
   aryuHH[15]=aryuHH[15]*MMt;
   aryuHH[16]=aryuHH[17]*aryuHH[14];
   aryuHH[16]=aryuHH[16] + aryuHH[11] - aryuHH[6];
   aryuHH[15]=aryuHH[15] + 4*aryuHH[16];
   aryuHH[15]=aryuHH[15]*aryuHH[4];
   aryuHH[16]= - 6 + 36*aryuHH[17];
   aryuHH[16]=aryuHH[16]*aryuHH[7];
   aryuHH[16]=aryuHH[16] + 5;
   aryuHH[17]=aryuHH[8]*aryuHH[4];
   aryuHH[17]=2*aryuHH[17];
   aryuHH[16]=aryuHH[16]*aryuHH[17];
   aryuHH[17]= - aryuHH[10]*pow(aryuHH[8],2)*aryuHH[4];
   aryuHH[15]=4*aryuHH[17] - aryuHH[15] - aryuHH[16];
   aryuHH[16]=pow(aryuHH[2],2);
   aryuHH[17]=pow(aryuHH[5],2);
   aryuHH[16]=aryuHH[16] + aryuHH[17];

      yuHHret = aryuHH[16]*aryuHH[15]*aryuHH[1];
      return yuHHret;
}
