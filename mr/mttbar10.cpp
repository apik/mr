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
tt::x10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbar[26], mttbarret;

    armttbar[1]=double(nH);
    armttbar[2]=double(boson);
    armttbar[3]=Tsil::A(MMt,mu2);
    armttbar[4]=pow(CW,-1);
    armttbar[5]=pow(MMH,-1);
    armttbar[6]=pow(MMZ,-1);
    armttbar[7]=pow(SW,-1);
    armttbar[8]=Tsil::A(MMb,mu2);
    armttbar[9]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbar[10]=Tsil::B(MMZ,MMt,MMt,mu2);
    armttbar[11]=pow(MMt,-1);
    armttbar[12]=Tsil::B(MMW,MMb,MMt,mu2);
    armttbar[13]=Tsil::A(MMH,mu2);
    armttbar[14]=Tsil::A(MMZ,mu2);
    armttbar[15]=Tsil::A(MMW,mu2);
   armttbar[16]=armttbar[8]*armttbar[11];
   armttbar[17]=1./2.*armttbar[16] - armttbar[12];
   armttbar[18]=1./8.*armttbar[12];
   armttbar[19]=armttbar[18]*MMb*armttbar[11];
   armttbar[20]=3*armttbar[5];
   armttbar[21]=armttbar[20]*armttbar[1];
   armttbar[22]= - armttbar[8]*armttbar[21];
   armttbar[23]=1./8.*armttbar[11];
   armttbar[23]= - armttbar[15]*armttbar[23];
   armttbar[17]=1./4.*armttbar[17] + armttbar[19] + armttbar[23] + 
   armttbar[22];
   armttbar[22]=pow(armttbar[4],2);
   armttbar[23]=pow(armttbar[7],2);
   armttbar[24]=armttbar[22] + armttbar[23];
   armttbar[17]=MMb*armttbar[24]*armttbar[17];
   armttbar[21]= - armttbar[3]*armttbar[21];
   armttbar[18]=armttbar[18] + armttbar[21];
   armttbar[18]=MMt*armttbar[18];
   armttbar[21]=armttbar[13] + 1./2.*armttbar[8] + armttbar[3];
   armttbar[25]= - 1./4.*MMH + MMt;
   armttbar[25]=armttbar[9]*armttbar[25];
   armttbar[18]=1./2.*armttbar[25] + 1./8.*armttbar[15] + 1./4.*
   armttbar[21] + armttbar[18];
   armttbar[18]=armttbar[24]*armttbar[18];
   armttbar[17]=armttbar[17] + armttbar[18];
   armttbar[17]=armttbar[6]*armttbar[17];
   armttbar[18]= - 1./2.*armttbar[11] + armttbar[20];
   armttbar[18]=armttbar[15]*armttbar[18];
   armttbar[16]=1./2.*armttbar[18] + armttbar[19] + 1./4.*armttbar[16]
    - 3./8.;
   armttbar[16]=armttbar[23]*armttbar[16];
   armttbar[18]=1./8.*armttbar[23];
   armttbar[19]=armttbar[18] + 17./72.*armttbar[22];
   armttbar[20]=armttbar[19] - 4./9.;
   armttbar[21]= - armttbar[14]*armttbar[20];
   armttbar[19]=8./9. + armttbar[19];
   armttbar[19]=armttbar[3]*armttbar[19];
   armttbar[19]=armttbar[19] + armttbar[21];
   armttbar[19]=armttbar[11]*armttbar[19];
   armttbar[21]=MMZ*armttbar[11];
   armttbar[20]= - armttbar[20]*armttbar[21];
   armttbar[18]=armttbar[20] - 7./72.*armttbar[22] + armttbar[18] + 8./
   9.;
   armttbar[18]=armttbar[10]*armttbar[18];
   armttbar[20]=armttbar[14]*armttbar[24];
   armttbar[24]=1./2.*armttbar[22] - 1 + 3./2.*armttbar[23];
   armttbar[24]=MMZ*armttbar[24];
   armttbar[20]=3./4.*armttbar[20] + armttbar[24];
   armttbar[20]=armttbar[5]*armttbar[20];
   armttbar[24]=1 - armttbar[23];
   armttbar[21]=armttbar[24]*armttbar[21];
   armttbar[21]=1./2.*armttbar[23] + armttbar[21];
   armttbar[21]=armttbar[12]*armttbar[21];
   armttbar[16]=armttbar[17] + armttbar[20] + armttbar[18] + 1./4.*
   armttbar[21] - 1./72.*armttbar[22] - 8./9. + armttbar[19] + 
   armttbar[16];

      mttbarret = armttbar[16]*armttbar[2];
      return mttbarret;
}
