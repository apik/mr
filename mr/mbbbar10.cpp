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
bb::x10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbar[24], mbbbarret;

    armbbbar[1]=double(nH);
    armbbbar[2]=double(boson);
    armbbbar[3]=Tsil::A(MMt,mu2);
    armbbbar[4]=pow(CW,-1);
    armbbbar[5]=pow(MMH,-1);
    armbbbar[6]=pow(MMZ,-1);
    armbbbar[7]=pow(SW,-1);
    armbbbar[8]=Tsil::A(MMb,mu2);
    armbbbar[9]=Tsil::B(MMH,MMb,MMb,mu2);
    armbbbar[10]=Tsil::B(MMZ,MMb,MMb,mu2);
    armbbbar[11]=pow(MMb,-1);
    armbbbar[12]=Tsil::B(MMW,MMt,MMb,mu2);
    armbbbar[13]=Tsil::A(MMH,mu2);
    armbbbar[14]=Tsil::A(MMZ,mu2);
    armbbbar[15]=Tsil::A(MMW,mu2);
   armbbbar[16]=3*armbbbar[5];
   armbbbar[17]=armbbbar[16]*armbbbar[1];
   armbbbar[18]=armbbbar[17]*MMt;
   armbbbar[19]=armbbbar[11]*MMt;
   armbbbar[19]=1./8.*armbbbar[19];
   armbbbar[18]= - armbbbar[18] + armbbbar[19] + 1./8.;
   armbbbar[18]=armbbbar[18]*armbbbar[3];
   armbbbar[20]= - MMb + 1./4.*MMH;
   armbbbar[21]=1./2.*armbbbar[9];
   armbbbar[20]=armbbbar[20]*armbbbar[21];
   armbbbar[21]=1./2.*armbbbar[15] + armbbbar[13] + armbbbar[8];
   armbbbar[19]=armbbbar[19]*armbbbar[15];
   armbbbar[17]=armbbbar[17]*armbbbar[8]*MMb;
   armbbbar[17]= - armbbbar[18] - 1./4.*armbbbar[21] + armbbbar[20] + 
   armbbbar[19] + armbbbar[17];
   armbbbar[17]=armbbbar[17]*armbbbar[6];
   armbbbar[18]= - 13 + 17*armbbbar[10];
   armbbbar[19]=MMZ + 3./2.*armbbbar[14];
   armbbbar[19]=armbbbar[5]*armbbbar[19];
   armbbbar[20]=MMZ*armbbbar[10];
   armbbbar[21]=armbbbar[20] + armbbbar[14];
   armbbbar[22]=armbbbar[8] - armbbbar[21];
   armbbbar[22]=armbbbar[11]*armbbbar[22];
   armbbbar[18]=5./36.*armbbbar[22] + 1./36.*armbbbar[18] + 
   armbbbar[19];
   armbbbar[18]= - armbbbar[17] + 1./2.*armbbbar[18];
   armbbbar[19]=pow(armbbbar[4],2);
   armbbbar[18]=armbbbar[19]*armbbbar[18];
   armbbbar[20]=armbbbar[20] - armbbbar[8];
   armbbbar[22]=armbbbar[15] + 1./2.*armbbbar[14];
   armbbbar[20]=armbbbar[3] - armbbbar[22] - 1./2.*armbbbar[20];
   armbbbar[23]=1./2.*armbbbar[11];
   armbbbar[20]=armbbbar[23]*armbbbar[20];
   armbbbar[22]=MMZ + armbbbar[22];
   armbbbar[16]=armbbbar[22]*armbbbar[16];
   armbbbar[22]= - 3 + armbbbar[10];
   armbbbar[16]=1./4.*armbbbar[22] + armbbbar[16] + armbbbar[20];
   armbbbar[20]=armbbbar[23]*pow(MMt,2);
   armbbbar[20]=1./2.*MMb + armbbbar[20] - MMt;
   armbbbar[20]=armbbbar[20]*armbbbar[6];
   armbbbar[22]= - MMZ + 1./2.*MMt;
   armbbbar[22]=armbbbar[11]*armbbbar[22];
   armbbbar[22]=armbbbar[20] + 1./2. + armbbbar[22];
   armbbbar[23]=1./4.*armbbbar[12];
   armbbbar[22]=armbbbar[22]*armbbbar[23];
   armbbbar[16]=armbbbar[22] + 1./2.*armbbbar[16] - armbbbar[17];
   armbbbar[16]=pow(armbbbar[7],2)*armbbbar[16];
   armbbbar[17]=armbbbar[19]*armbbbar[20];
   armbbbar[19]=armbbbar[11]*MMZ;
   armbbbar[17]=armbbbar[19] + armbbbar[17];
   armbbbar[17]=armbbbar[17]*armbbbar[23];
   armbbbar[19]= - 1 + armbbbar[10];
   armbbbar[20]= - armbbbar[5]*MMZ;
   armbbbar[21]=2*armbbbar[8] + armbbbar[21];
   armbbbar[21]=armbbbar[11]*armbbbar[21];
   armbbbar[16]=armbbbar[16] + armbbbar[17] + 1./9.*armbbbar[21] + 2./9.
   *armbbbar[19] + armbbbar[20] + armbbbar[18];

      mbbbarret = armbbbar[16]*armbbbar[2];
      return mbbbarret;
}
