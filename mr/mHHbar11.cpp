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
HH<OS>::x11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbar[20], mHHbarret;

    armHHbar[1]=double(nH);
    armHHbar[2]=pow(CW,-1);
    armHHbar[3]=pow(MMH,-1);
    armHHbar[4]=pow(MMZ,-1);
    armHHbar[5]=pow(SW,-1);
    armHHbar[6]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbar[7]=Tsil::A(MMt,mu2);
    armHHbar[8]=Tsil::Beps(MMt,MMt,MMH,mu2);
    armHHbar[9]=Tsil::Aeps(MMt,mu2);
    armHHbar[10]=prottttt0->M(0);
    armHHbar[11]=prottttt0->Vxzuv(0);
    armHHbar[12]=prottttt0->Suxv(0);
   armHHbar[13]=2*armHHbar[11];
   armHHbar[14]=armHHbar[13] + armHHbar[10];
   armHHbar[15]=8*MMt;
   armHHbar[14]=armHHbar[14]*armHHbar[15];
   armHHbar[14]=4*armHHbar[8] + armHHbar[14] - 9;
   armHHbar[15]=pow(MMt,2);
   armHHbar[14]=armHHbar[15]*armHHbar[14];
   armHHbar[16]=6*MMt;
   armHHbar[17]=armHHbar[16]*armHHbar[9];
   armHHbar[18]=armHHbar[12]*MMt;
   armHHbar[19]=pow(armHHbar[7],2);
   armHHbar[14]= - armHHbar[17] - armHHbar[14] + armHHbar[18] + 12*
   armHHbar[19];
   armHHbar[17]=armHHbar[3]*armHHbar[1];
   armHHbar[14]=armHHbar[17]*armHHbar[14];
   armHHbar[13]=armHHbar[13] + 3*armHHbar[10];
   armHHbar[13]=armHHbar[13]*MMt;
   armHHbar[13]=armHHbar[13] - 3;
   armHHbar[18]=MMH*armHHbar[10];
   armHHbar[13]= - armHHbar[18] + 2*armHHbar[13];
   armHHbar[18]=armHHbar[1]*MMt;
   armHHbar[13]=armHHbar[13]*armHHbar[18];
   armHHbar[13]=armHHbar[13] + armHHbar[14];
   armHHbar[14]=armHHbar[15]*armHHbar[17];
   armHHbar[15]= - armHHbar[18] + 5*armHHbar[14];
   armHHbar[16]=armHHbar[16]*armHHbar[17];
   armHHbar[16]=armHHbar[16] - armHHbar[1];
   armHHbar[17]=3*armHHbar[7];
   armHHbar[16]=armHHbar[16]*armHHbar[17];
   armHHbar[15]=armHHbar[16] + 2*armHHbar[15];
   armHHbar[14]= - 3*armHHbar[18] + 14*armHHbar[14];
   armHHbar[14]=armHHbar[14]*armHHbar[6];
   armHHbar[14]=armHHbar[14] + 2*armHHbar[15];
   armHHbar[14]=armHHbar[14]*armHHbar[6];
   armHHbar[13]=armHHbar[14] + 2*armHHbar[13];
   armHHbar[14]= - pow(armHHbar[5],2);
   armHHbar[15]= - pow(armHHbar[2],2);
   armHHbar[14]=armHHbar[14] + armHHbar[15];

      mHHbarret = 2*armHHbar[14]*armHHbar[13]*armHHbar[4];
      return mHHbarret;
}
