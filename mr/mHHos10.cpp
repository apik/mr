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
HH<MS>::x10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHos[22], mHHosret;

    armHHos[1]=double(nH);
    armHHos[2]=Tsil::B(mmt,mmt,mmH,mu2);
    armHHos[3]=pow(mmZ,-1);
    armHHos[4]=pow(mmH,-1);
    armHHos[5]=pow(s,-1);
    armHHos[6]=pow(c,-1);
    armHHos[7]=Tsil::A(mmt,mu2);
    armHHos[8]=double(boson);
    armHHos[9]=Tsil::B(mmW,mmW,mmH,mu2);
    armHHos[10]=Tsil::B(mmZ,mmZ,mmH,mu2);
    armHHos[11]=Tsil::B(mmH,mmH,mmH,mu2);
    armHHos[12]=Tsil::A(mmW,mu2);
    armHHos[13]=Tsil::A(mmZ,mu2);
    armHHos[14]=Tsil::A(mmH,mu2);
   armHHos[15]=pow(armHHos[6],2);
   armHHos[16]=pow(armHHos[5],2);
   armHHos[15]=armHHos[15] + armHHos[16];
   armHHos[17]= - 9./2.*armHHos[11] - armHHos[9];
   armHHos[17]=armHHos[15]*armHHos[17];
   armHHos[18]=1./2.*armHHos[10];
   armHHos[18]=armHHos[15]*armHHos[18];
   armHHos[17]= - armHHos[18] + armHHos[17];
   armHHos[17]=mmH*armHHos[17];
   armHHos[19]=armHHos[15]*armHHos[13];
   armHHos[17]= - armHHos[19] + armHHos[17];
   armHHos[20]=armHHos[14]*armHHos[15];
   armHHos[17]= - 3./4.*armHHos[20] + 1./4.*armHHos[17];
   armHHos[17]=armHHos[3]*armHHos[17];
   armHHos[20]=1 - armHHos[16];
   armHHos[20]=armHHos[9]*armHHos[20];
   armHHos[20]= - armHHos[18] + armHHos[20];
   armHHos[20]=mmZ*armHHos[20];
   armHHos[21]= - armHHos[12]*armHHos[16];
   armHHos[19]=armHHos[21] - 1./2.*armHHos[19] + armHHos[20];
   armHHos[19]=armHHos[4]*armHHos[19];
   armHHos[20]=1./2.*armHHos[3];
   armHHos[20]= - armHHos[12]*armHHos[15]*armHHos[20];
   armHHos[16]=armHHos[9]*armHHos[16];
   armHHos[16]=3*armHHos[19] + armHHos[20] + armHHos[16] + armHHos[18]
    + armHHos[17];
   armHHos[16]=armHHos[8]*armHHos[16];
   armHHos[17]=armHHos[15]*mmt*armHHos[2];
   armHHos[15]=armHHos[7]*armHHos[15];
   armHHos[15]=armHHos[15] + armHHos[17];
   armHHos[15]=armHHos[4]*mmt*armHHos[15];
   armHHos[15]= - 1./2.*armHHos[17] + 2*armHHos[15];
   armHHos[15]=armHHos[15]*armHHos[3]*armHHos[1];

      mHHosret = 3*armHHos[15] + armHHos[16];
      return mHHosret;
}
