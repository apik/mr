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
bb::my11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubb[36], yubbret;

    aryubb[1]=double(boson);
    aryubb[2]=pow(CW,-1);
    aryubb[3]=pow(MMZ,-1);
    aryubb[4]=pow(SW,-1);
    aryubb[5]=Tsil::I2(0,MMW,MMt,mu2);
    aryubb[6]=Tsil::I2(0,0,MMt,mu2);
    aryubb[7]=Tsil::A(MMZ,mu2);
    aryubb[8]=Tsil::A(MMb,mu2);
    aryubb[9]=pow(MMb,-1);
    aryubb[10]=Tsil::A(MMW,mu2);
    aryubb[11]=Tsil::A(MMt,mu2);
    aryubb[12]=pow(MMt,-1);
    aryubb[13]=Tsil::Aeps(MMW,mu2);
    aryubb[14]=Tsil::Aeps(MMt,mu2);
    aryubb[15]=Tsil::Aeps(MMb,mu2);
    aryubb[16]=prot0bb0b->Tvxu(0);
    aryubb[17]=1/(MMt - MMW);
    aryubb[18]=1/( - MMW + MMH);
    aryubb[19]=Tsil::A(MMH,mu2);
   aryubb[20]=aryubb[8]*aryubb[9];
   aryubb[21]=3*aryubb[20];
   aryubb[22]=aryubb[21] - 1;
   aryubb[23]=aryubb[7] - aryubb[10];
   aryubb[24]=pow(aryubb[4],2);
   aryubb[25]=1./2.*aryubb[24];
   aryubb[23]= - aryubb[25]*aryubb[23]*aryubb[22];
   aryubb[26]= - aryubb[7] + aryubb[14] + aryubb[6];
   aryubb[27]=aryubb[13] - aryubb[5];
   aryubb[26]= - 1./12.*MMH - aryubb[10] + 7./6.*MMt + 4*aryubb[27] + 2
   *aryubb[26];
   aryubb[28]=3*aryubb[10];
   aryubb[29]=aryubb[28] - 5./4.*MMt + 1./4.*MMH;
   aryubb[30]=3./2.*aryubb[7] + aryubb[29];
   aryubb[30]=aryubb[30]*aryubb[20];
   aryubb[21]=aryubb[21] - 5;
   aryubb[31]=aryubb[11]*aryubb[12];
   aryubb[31]= - 1./2.*aryubb[21] - 2*aryubb[31];
   aryubb[31]=aryubb[11]*aryubb[31];
   aryubb[23]=aryubb[23] + aryubb[31] + aryubb[30] + aryubb[26];
   aryubb[23]=aryubb[23]*aryubb[24];
   aryubb[30]=pow(aryubb[2],2);
   aryubb[26]=aryubb[30]*aryubb[26];
   aryubb[31]=aryubb[30]*aryubb[9];
   aryubb[29]=5./6.*aryubb[7] + aryubb[29];
   aryubb[29]=aryubb[29]*aryubb[31];
   aryubb[32]=aryubb[9]*aryubb[7];
   aryubb[29]= - 4./3.*aryubb[32] + aryubb[29];
   aryubb[29]=aryubb[8]*aryubb[29];
   aryubb[21]= - aryubb[30]*aryubb[21];
   aryubb[32]=aryubb[30]*aryubb[11];
   aryubb[33]=aryubb[12]*aryubb[32];
   aryubb[21]=1./2.*aryubb[21] - 2*aryubb[33];
   aryubb[21]=aryubb[11]*aryubb[21];
   aryubb[21]=aryubb[23] + aryubb[21] + aryubb[29] + aryubb[26];
   aryubb[21]=aryubb[3]*aryubb[21];
   aryubb[23]=3./2.*aryubb[20];
   aryubb[26]=aryubb[23]*aryubb[10];
   aryubb[26]=aryubb[26] + aryubb[28];
   aryubb[27]=aryubb[27] + aryubb[14];
   aryubb[23]=aryubb[23] - 11;
   aryubb[23]=aryubb[23]*aryubb[11];
   aryubb[23]= - aryubb[23] + aryubb[26] + 7*aryubb[27];
   aryubb[28]=aryubb[24] - 1;
   aryubb[23]=aryubb[23]*aryubb[28];
   aryubb[29]=pow(CW,2);
   aryubb[33]=aryubb[29] + 1;
   aryubb[34]=aryubb[33] - aryubb[24];
   aryubb[35]=aryubb[11] + aryubb[10] + aryubb[27];
   aryubb[35]= - aryubb[35]*aryubb[34];
   aryubb[29]=aryubb[33]*aryubb[29];
   aryubb[29]=aryubb[29] - aryubb[28];
   aryubb[29]=MMZ*aryubb[29];
   aryubb[29]=2*aryubb[29] + aryubb[35];
   aryubb[29]=MMZ*aryubb[29];
   aryubb[33]=aryubb[24]*aryubb[11];
   aryubb[35]=aryubb[11] - aryubb[33];
   aryubb[35]=aryubb[10]*aryubb[35];
   aryubb[29]=aryubb[29] + aryubb[35];
   aryubb[29]=aryubb[17]*aryubb[29];
   aryubb[34]=MMZ*aryubb[34];
   aryubb[23]=3*aryubb[29] + 17*aryubb[34] + aryubb[23];
   aryubb[23]=MMZ*aryubb[23];
   aryubb[29]=4*aryubb[10];
   aryubb[34]= - aryubb[29] - 3./2.*aryubb[11];
   aryubb[34]=aryubb[34]*aryubb[33];
   aryubb[23]=aryubb[34] + aryubb[23];
   aryubb[23]=aryubb[17]*aryubb[23];
   aryubb[34]= - aryubb[30]*aryubb[29];
   aryubb[32]=aryubb[34] - 1./2.*aryubb[32];
   aryubb[32]=aryubb[11]*aryubb[32];
   aryubb[29]= - aryubb[29] - 1./2.*aryubb[11];
   aryubb[29]=aryubb[29]*aryubb[33];
   aryubb[29]=aryubb[32] + aryubb[29];
   aryubb[29]=aryubb[3]*aryubb[29];
   aryubb[26]=9*aryubb[11] + aryubb[26] + 8*aryubb[27];
   aryubb[24]=aryubb[26]*aryubb[24];
   aryubb[20]= - 10 + 1./2.*aryubb[20];
   aryubb[20]=MMZ*aryubb[20]*aryubb[28];
   aryubb[20]=aryubb[23] + aryubb[29] + aryubb[24] + 3*aryubb[20];
   aryubb[20]=aryubb[17]*aryubb[20];
   aryubb[23]=aryubb[8]*pow(aryubb[9],2);
   aryubb[23]=32*aryubb[23] + 34*aryubb[9] - 11./2.*aryubb[31];
   aryubb[23]=aryubb[8]*aryubb[23];
   aryubb[24]=aryubb[9]*aryubb[15];
   aryubb[24]=aryubb[24] + aryubb[16];
   aryubb[24]= - 77./3. - 20*aryubb[24];
   aryubb[23]=aryubb[23] + 2*aryubb[24] - 125./24.*aryubb[30];
   aryubb[24]=aryubb[10] - aryubb[19];
   aryubb[22]= - aryubb[18]*aryubb[24]*aryubb[22];
   aryubb[22]= - 169./4. + aryubb[22];
   aryubb[22]=aryubb[22]*aryubb[25];
   aryubb[20]=aryubb[20] + aryubb[21] + 1./9.*aryubb[23] + aryubb[22];

      yubbret = aryubb[20]*aryubb[1];
      return yubbret;
}
