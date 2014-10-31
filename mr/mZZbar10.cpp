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

#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZbar[31], mZZbarret;

    armZZbar[1]=double(nH);
    armZZbar[2]=double(boson);
    armZZbar[3]=pow(CW,-1);
    armZZbar[4]=pow(MMZ,-1);
    armZZbar[5]=pow(SW,-1);
    armZZbar[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    armZZbar[7]=Tsil::B(MMb,MMb,MMZ,mu2);
    armZZbar[8]=Tsil::B(0,0,MMZ,mu2);
    armZZbar[9]=Tsil::A(MMt,mu2);
    armZZbar[10]=pow(MMH,-1);
    armZZbar[11]=Tsil::A(MMb,mu2);
    armZZbar[12]=double(nL + nH);
    armZZbar[13]=Tsil::B(MMZ,MMH,MMZ,mu2);
    armZZbar[14]=Tsil::B(MMW,MMW,MMZ,mu2);
    armZZbar[15]=Tsil::A(MMH,mu2);
    armZZbar[16]=Tsil::A(MMZ,mu2);
    armZZbar[17]=Tsil::A(MMW,mu2);
   armZZbar[18]=pow(armZZbar[3],2);
   armZZbar[19]=pow(armZZbar[5],2);
   armZZbar[20]=5./9.*armZZbar[18] + armZZbar[19] - 8./9.;
   armZZbar[21]=1./2.*armZZbar[19];
   armZZbar[22]=armZZbar[21] + 17./18.*armZZbar[18];
   armZZbar[23]= - 8./9. - armZZbar[22];
   armZZbar[23]=armZZbar[7]*armZZbar[23];
   armZZbar[24]=armZZbar[18] + armZZbar[19];
   armZZbar[25]=6*armZZbar[10];
   armZZbar[25]=armZZbar[24]*armZZbar[25];
   armZZbar[26]= - armZZbar[11]*armZZbar[25];
   armZZbar[23]=armZZbar[26] + armZZbar[23] + armZZbar[20];
   armZZbar[23]=MMb*armZZbar[23];
   armZZbar[26]=armZZbar[19] - 32./9. + 17./9.*armZZbar[18];
   armZZbar[27]=7./18.*armZZbar[18] - 32./9. - armZZbar[21];
   armZZbar[27]=armZZbar[6]*armZZbar[27];
   armZZbar[27]=armZZbar[27] + armZZbar[26];
   armZZbar[27]=MMt*armZZbar[27];
   armZZbar[20]=armZZbar[11]*armZZbar[20];
   armZZbar[20]=armZZbar[23] + armZZbar[20] + armZZbar[27];
   armZZbar[20]=armZZbar[4]*armZZbar[20];
   armZZbar[23]= - MMt*armZZbar[25];
   armZZbar[23]=armZZbar[23] + armZZbar[26];
   armZZbar[23]=armZZbar[9]*armZZbar[4]*armZZbar[23];
   armZZbar[21]=5./18.*armZZbar[18] - 4./9. + armZZbar[21];
   armZZbar[21]=armZZbar[7]*armZZbar[21];
   armZZbar[22]= - 16./9. + armZZbar[22];
   armZZbar[22]=armZZbar[6]*armZZbar[22];
   armZZbar[25]= - 11./9.*armZZbar[18] + 20./9. - armZZbar[19];
   armZZbar[25]=armZZbar[8]*armZZbar[25];
   armZZbar[20]=armZZbar[23] + armZZbar[20] + armZZbar[25] + 
   armZZbar[21] + armZZbar[22];
   armZZbar[20]=armZZbar[1]*armZZbar[20];
   armZZbar[21]=armZZbar[16]*armZZbar[24];
   armZZbar[22]=3*armZZbar[19];
   armZZbar[23]=armZZbar[18] - 2 + armZZbar[22];
   armZZbar[23]=MMZ*armZZbar[23];
   armZZbar[21]=armZZbar[23] + 3./2.*armZZbar[21];
   armZZbar[21]=armZZbar[10]*armZZbar[21];
   armZZbar[23]=armZZbar[13] + 1./2.;
   armZZbar[23]= - armZZbar[23]*armZZbar[24];
   armZZbar[25]=armZZbar[24]*armZZbar[4];
   armZZbar[26]=armZZbar[15] - armZZbar[16];
   armZZbar[26]=armZZbar[26]*armZZbar[25];
   armZZbar[23]=1./4.*armZZbar[26] + armZZbar[23];
   armZZbar[23]=armZZbar[4]*armZZbar[23];
   armZZbar[24]=MMH*pow(armZZbar[4],2)*armZZbar[13]*armZZbar[24];
   armZZbar[23]=armZZbar[23] + 1./4.*armZZbar[24];
   armZZbar[23]=MMH*armZZbar[23];
   armZZbar[24]=pow(CW,2);
   armZZbar[24]=4*armZZbar[24];
   armZZbar[26]=1./12.*armZZbar[18] + armZZbar[24] + 29./3. - 33./4.*
   armZZbar[19];
   armZZbar[26]=armZZbar[14]*armZZbar[26];
   armZZbar[27]= - 59./18. + armZZbar[13];
   armZZbar[27]=armZZbar[27]*armZZbar[19];
   armZZbar[28]= - 1./18. + armZZbar[13];
   armZZbar[28]=armZZbar[28]*armZZbar[18];
   armZZbar[29]=armZZbar[15] + 1./3.*armZZbar[16];
   armZZbar[25]=armZZbar[29]*armZZbar[25];
   armZZbar[29]=5./3.*armZZbar[18] + armZZbar[19] - 8./3.;
   armZZbar[30]= - 1./3. + armZZbar[8];
   armZZbar[29]=armZZbar[12]*armZZbar[29]*armZZbar[30];
   armZZbar[18]=1./6.*armZZbar[18] + 4 - 5./2.*armZZbar[19];
   armZZbar[18]=armZZbar[4]*armZZbar[18];
   armZZbar[19]=armZZbar[10]*armZZbar[22];
   armZZbar[18]=armZZbar[19] + armZZbar[18];
   armZZbar[18]=armZZbar[17]*armZZbar[18];
   armZZbar[18]=armZZbar[26] + 1./3.*armZZbar[23] + armZZbar[18] + 4./3.
   *armZZbar[29] + 1./2.*armZZbar[25] + armZZbar[28] + armZZbar[24] + 8.
   /3. + armZZbar[27] + armZZbar[21] + armZZbar[20];

      mZZbarret = armZZbar[18]*armZZbar[2];
      return mZZbarret;
}
