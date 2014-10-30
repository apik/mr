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
HH<OS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbar[24], mHHbarret;

    armHHbar[1]=double(nH);
    armHHbar[2]=double(boson);
    armHHbar[3]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbar[4]=pow(CW,-1);
    armHHbar[5]=pow(MMH,-1);
    armHHbar[6]=pow(MMZ,-1);
    armHHbar[7]=pow(SW,-1);
    armHHbar[8]=Tsil::B(MMb,MMb,MMH,mu2);
    armHHbar[9]=Tsil::A(MMt,mu2);
    armHHbar[10]=Tsil::A(MMb,mu2);
    armHHbar[11]=Tsil::B(MMH,MMH,MMH,mu2);
    armHHbar[12]=Tsil::B(MMZ,MMZ,MMH,mu2);
    armHHbar[13]=Tsil::B(MMW,MMW,MMH,mu2);
    armHHbar[14]=Tsil::A(MMH,mu2);
    armHHbar[15]=Tsil::A(MMZ,mu2);
    armHHbar[16]=Tsil::A(MMW,mu2);
   armHHbar[17]= - MMt*armHHbar[3];
   armHHbar[17]=armHHbar[17] - armHHbar[9];
   armHHbar[18]=2*armHHbar[5];
   armHHbar[17]=armHHbar[18]*armHHbar[17];
   armHHbar[17]=1./2.*armHHbar[3] + armHHbar[17];
   armHHbar[19]=pow(armHHbar[4],2);
   armHHbar[20]=pow(armHHbar[7],2);
   armHHbar[19]=armHHbar[19] + armHHbar[20];
   armHHbar[21]=armHHbar[19]*armHHbar[1];
   armHHbar[17]=MMt*armHHbar[21]*armHHbar[17];
   armHHbar[21]=armHHbar[21]*MMb;
   armHHbar[18]= - MMb*armHHbar[18];
   armHHbar[18]=1./2. + armHHbar[18];
   armHHbar[18]=armHHbar[8]*armHHbar[18]*armHHbar[21];
   armHHbar[17]=armHHbar[17] + armHHbar[18];
   armHHbar[18]=armHHbar[19]*armHHbar[15];
   armHHbar[18]=1./2.*armHHbar[18];
   armHHbar[22]=armHHbar[16]*armHHbar[19];
   armHHbar[22]=armHHbar[18] + armHHbar[22];
   armHHbar[23]=1./4.*armHHbar[13] + 1./8.*armHHbar[12] + 9./8.*
   armHHbar[11];
   armHHbar[23]=armHHbar[23]*MMH;
   armHHbar[23]=3./4.*armHHbar[14] + armHHbar[23];
   armHHbar[23]=armHHbar[19]*armHHbar[23];
   armHHbar[21]=armHHbar[10]*armHHbar[5]*armHHbar[21];
   armHHbar[17]= - 6*armHHbar[21] + 1./2.*armHHbar[22] + 3*armHHbar[17]
    + armHHbar[23];
   armHHbar[17]=armHHbar[6]*armHHbar[17];
   armHHbar[21]=armHHbar[16]*armHHbar[20];
   armHHbar[18]=armHHbar[18] + armHHbar[21];
   armHHbar[21]=3*armHHbar[5];
   armHHbar[18]=armHHbar[18]*armHHbar[21];
   armHHbar[21]=armHHbar[21]*MMZ;
   armHHbar[22]=armHHbar[21] - 1;
   armHHbar[19]=armHHbar[12]*armHHbar[19]*armHHbar[22];
   armHHbar[22]= - 1 + armHHbar[20];
   armHHbar[21]=armHHbar[22]*armHHbar[21];
   armHHbar[20]= - armHHbar[20] + armHHbar[21];
   armHHbar[20]=armHHbar[13]*armHHbar[20];
   armHHbar[17]=armHHbar[17] + armHHbar[20] + armHHbar[18] + 1./2.*
   armHHbar[19];

      mHHbarret = armHHbar[17]*armHHbar[2];
      return mHHbarret;
}
