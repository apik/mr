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
HH<MS>::x11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHos[19], mHHosret;

    armHHos[1]=double(nH);
    armHHos[2]=pow(mmZ,-1);
    armHHos[3]=pow(mmH,-1);
    armHHos[4]=pow(s,-1);
    armHHos[5]=pow(c,-1);
    armHHos[6]=Tsil::B(mmt,mmt,mmH,mu2);
    armHHos[7]=Tsil::A(mmt,mu2);
    armHHos[8]=Tsil::Beps(mmt,mmt,mmH,mu2);
    armHHos[9]=Tsil::Aeps(mmt,mu2);
    armHHos[10]=prottttt0->M(0);
    armHHos[11]=prottttt0->Vxzuv(0);
    armHHos[12]=prottttt0->Suxv(0);
   armHHos[13]=2*armHHos[11];
   armHHos[14]=armHHos[13] + armHHos[10];
   armHHos[15]=8*mmt;
   armHHos[14]=armHHos[14]*armHHos[15];
   armHHos[15]=20 + 7*armHHos[6];
   armHHos[15]=armHHos[15]*armHHos[6];
   armHHos[14]= - armHHos[14] + armHHos[15] + 11 - 4*armHHos[8];
   armHHos[15]=armHHos[3]*armHHos[1];
   armHHos[14]=armHHos[14]*armHHos[15];
   armHHos[13]=armHHos[13] + 3*armHHos[10];
   armHHos[16]=2*armHHos[1];
   armHHos[13]=armHHos[13]*armHHos[16];
   armHHos[13]=armHHos[14] + armHHos[13];
   armHHos[14]=2*mmt;
   armHHos[13]=armHHos[13]*armHHos[14];
   armHHos[14]= - 1 + 3*armHHos[6];
   armHHos[17]=4*armHHos[7];
   armHHos[14]=armHHos[14]*armHHos[17];
   armHHos[14]=armHHos[14] - armHHos[12] + 6*armHHos[9];
   armHHos[17]=2*armHHos[15];
   armHHos[14]=armHHos[14]*armHHos[17];
   armHHos[17]=armHHos[6] + 2;
   armHHos[17]=armHHos[17]*armHHos[6];
   armHHos[17]=armHHos[17] + 4;
   armHHos[18]=3*armHHos[1];
   armHHos[17]=armHHos[17]*armHHos[18];
   armHHos[16]=armHHos[10]*armHHos[16]*mmH;
   armHHos[13]= - armHHos[16] + armHHos[13] - armHHos[14] - armHHos[17]
   ;
   armHHos[13]=mmt*armHHos[13];
   armHHos[14]=armHHos[15]*pow(armHHos[7],2);
   armHHos[13]= - 36*armHHos[14] + armHHos[13];
   armHHos[14]=pow(armHHos[4],2);
   armHHos[15]=pow(armHHos[5],2);
   armHHos[14]=armHHos[14] + armHHos[15];

      mHHosret = 2*armHHos[14]*armHHos[13]*armHHos[2];
      return mHHosret;
}
