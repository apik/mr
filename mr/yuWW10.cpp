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

#include <WW.hpp>
std::complex<long double>
WW<OS>::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuWW[28], yuWWret;

    aryuWW[1]=double(nL + nH);
    aryuWW[2]=pow(SW,-1);
    aryuWW[3]=std::real(Tsil::B(0,0,MMW,mu2));
    aryuWW[4]=double(nH);
    aryuWW[5]=pow(CW,-1);
    aryuWW[6]=pow(MMZ,-1);
    aryuWW[7]=Tsil::B(0,MMt,MMW,mu2);
    aryuWW[8]=Tsil::A(MMt,mu2);
    aryuWW[9]=double(nL);
    aryuWW[10]=double(boson);
    aryuWW[11]=Tsil::B(MMW,MMH,MMW,mu2);
    aryuWW[12]=Tsil::B(MMW,MMZ,MMW,mu2);
    aryuWW[13]=Tsil::A(MMH,mu2);
    aryuWW[14]=Tsil::A(MMZ,mu2);
    aryuWW[15]=Tsil::A(MMW,mu2);
    aryuWW[16]=1/( - MMW + MMH);
   aryuWW[17]=MMH*aryuWW[10];
   aryuWW[18]=aryuWW[17]*aryuWW[11];
   aryuWW[19]=aryuWW[15] - aryuWW[13];
   aryuWW[20]=aryuWW[19]*aryuWW[10];
   aryuWW[18]=aryuWW[18] - aryuWW[20];
   aryuWW[20]=1./6.*MMH;
   aryuWW[18]=aryuWW[18]*aryuWW[20];
   aryuWW[20]=MMt*aryuWW[4];
   aryuWW[21]=aryuWW[20]*aryuWW[7];
   aryuWW[22]=aryuWW[4]*aryuWW[8];
   aryuWW[21]=aryuWW[21] + aryuWW[22];
   aryuWW[21]=aryuWW[21]*MMt;
   aryuWW[18]=aryuWW[18] - aryuWW[21];
   aryuWW[21]=1./2.*aryuWW[6];
   aryuWW[21]=aryuWW[18]*aryuWW[21];
   aryuWW[23]=aryuWW[11] + 1./8.;
   aryuWW[17]=1./3.*aryuWW[17];
   aryuWW[17]=aryuWW[23]*aryuWW[17];
   aryuWW[17]=aryuWW[21] - aryuWW[17];
   aryuWW[21]=aryuWW[7] - 1./2.;
   aryuWW[20]=aryuWW[21]*aryuWW[20];
   aryuWW[20]=aryuWW[20] + aryuWW[22];
   aryuWW[21]=1./2.*aryuWW[13];
   aryuWW[22]=5./6.*aryuWW[15] - aryuWW[14] - aryuWW[21];
   aryuWW[22]=aryuWW[10]*aryuWW[22];
   aryuWW[22]=aryuWW[22] - aryuWW[20];
   aryuWW[22]=1./2.*aryuWW[22] + aryuWW[17];
   aryuWW[22]=aryuWW[6]*aryuWW[22];
   aryuWW[23]=aryuWW[15] - aryuWW[14];
   aryuWW[24]=pow(aryuWW[2],2);
   aryuWW[25]=aryuWW[6]*aryuWW[23]*aryuWW[24];
   aryuWW[19]=aryuWW[16]*aryuWW[19];
   aryuWW[19]=aryuWW[19] - aryuWW[25];
   aryuWW[19]= - 209./72. + aryuWW[11] - 3./4.*aryuWW[19];
   aryuWW[19]=aryuWW[10]*aryuWW[19];
   aryuWW[25]=aryuWW[12]*aryuWW[10];
   aryuWW[26]= - 1./3. + aryuWW[7];
   aryuWW[26]=aryuWW[4]*aryuWW[26];
   aryuWW[27]=aryuWW[9] + 1./3.*aryuWW[1];
   aryuWW[27]=aryuWW[3]*aryuWW[27];
   aryuWW[19]= - 33./4.*aryuWW[25] + aryuWW[27] - 1./9.*aryuWW[1] + 
   aryuWW[22] - 1./3.*aryuWW[9] + aryuWW[26] + aryuWW[19];
   aryuWW[19]=aryuWW[19]*aryuWW[24];
   aryuWW[21]=53./6.*aryuWW[15] + 3*aryuWW[14] - aryuWW[21];
   aryuWW[21]=aryuWW[10]*aryuWW[21];
   aryuWW[20]=aryuWW[21] - aryuWW[20];
   aryuWW[17]=1./2.*aryuWW[20] + aryuWW[17];
   aryuWW[17]=aryuWW[6]*aryuWW[17];
   aryuWW[18]=aryuWW[6]*aryuWW[18];
   aryuWW[20]=aryuWW[10]*aryuWW[23];
   aryuWW[18]= - 1./6.*aryuWW[20] + aryuWW[18];
   aryuWW[18]=aryuWW[6]*aryuWW[18];
   aryuWW[18]=aryuWW[18] + 1./6.*aryuWW[25];
   aryuWW[20]=pow(aryuWW[5],2);
   aryuWW[18]=aryuWW[18]*aryuWW[20];
   aryuWW[17]=1./2.*aryuWW[18] + 17./12.*aryuWW[25] - 1./24.*aryuWW[10]
    + aryuWW[17];
   aryuWW[17]=aryuWW[17]*aryuWW[20];
   aryuWW[18]= - aryuWW[10] + aryuWW[25];

      yuWWret = aryuWW[17] + 4*aryuWW[18] + aryuWW[19];
      return yuWWret;
}
