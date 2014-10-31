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
WW<OS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWbar[29], mWWbarret;

    armWWbar[1]=double(nH);
    armWWbar[2]=double(boson);
    armWWbar[3]=pow(CW,-1);
    armWWbar[4]=pow(MMZ,-1);
    armWWbar[5]=pow(SW,-1);
    armWWbar[6]=Tsil::B(MMt,MMb,MMW,mu2);
    armWWbar[7]=Tsil::B(0,0,MMW,mu2);
    armWWbar[8]=Tsil::A(MMt,mu2);
    armWWbar[9]=pow(MMH,-1);
    armWWbar[10]=Tsil::A(MMb,mu2);
    armWWbar[11]=double(nL + nH);
    armWWbar[12]=Tsil::B(MMW,MMH,MMW,mu2);
    armWWbar[13]=Tsil::B(MMW,MMZ,MMW,mu2);
    armWWbar[14]=Tsil::A(MMH,mu2);
    armWWbar[15]=Tsil::A(MMZ,mu2);
    armWWbar[16]=Tsil::A(MMW,mu2);
   armWWbar[17]=armWWbar[9]*armWWbar[1];
   armWWbar[17]=6*armWWbar[17];
   armWWbar[18]=armWWbar[17]*MMb;
   armWWbar[18]=armWWbar[18] - armWWbar[1];
   armWWbar[19]=armWWbar[10]*armWWbar[2];
   armWWbar[18]=armWWbar[18]*armWWbar[19];
   armWWbar[20]=armWWbar[12]*armWWbar[2];
   armWWbar[21]=armWWbar[20] + 1./2.*armWWbar[2];
   armWWbar[22]=1./3.*MMH;
   armWWbar[21]=armWWbar[21]*armWWbar[22];
   armWWbar[18]=armWWbar[18] + armWWbar[21];
   armWWbar[21]=1./2.*armWWbar[1];
   armWWbar[22]=armWWbar[6]*armWWbar[21];
   armWWbar[22]=armWWbar[1] - armWWbar[22];
   armWWbar[23]=MMb + MMt;
   armWWbar[22]=armWWbar[23]*armWWbar[22];
   armWWbar[17]=armWWbar[17]*MMt;
   armWWbar[17]=armWWbar[17] - armWWbar[1];
   armWWbar[17]=armWWbar[17]*armWWbar[8];
   armWWbar[17]= - armWWbar[17] + 1./2.*armWWbar[14] + armWWbar[22];
   armWWbar[22]=pow(armWWbar[3],2);
   armWWbar[23]=1./12.*armWWbar[22];
   armWWbar[24]= - armWWbar[16] + armWWbar[15];
   armWWbar[24]=armWWbar[24]*armWWbar[23];
   armWWbar[24]=armWWbar[24] + 3./4.*armWWbar[15] + 35./12.*
   armWWbar[16] + armWWbar[17];
   armWWbar[24]=armWWbar[2]*armWWbar[24];
   armWWbar[24]=armWWbar[24] - armWWbar[18];
   armWWbar[24]=armWWbar[24]*armWWbar[22];
   armWWbar[25]= - MMt + 1./2.*MMb;
   armWWbar[25]=armWWbar[25]*MMb;
   armWWbar[26]=pow(MMt,2);
   armWWbar[25]=armWWbar[25] + 1./2.*armWWbar[26];
   armWWbar[26]=armWWbar[6]*armWWbar[1];
   armWWbar[25]=armWWbar[25]*armWWbar[26];
   armWWbar[27]=MMb - MMt;
   armWWbar[21]=armWWbar[27]*armWWbar[21];
   armWWbar[27]=armWWbar[21]*armWWbar[8];
   armWWbar[25]=armWWbar[25] - armWWbar[27];
   armWWbar[25]=armWWbar[25]*armWWbar[2];
   armWWbar[27]=armWWbar[14] - armWWbar[16];
   armWWbar[27]=armWWbar[27]*armWWbar[2];
   armWWbar[28]=armWWbar[20]*MMH;
   armWWbar[27]=armWWbar[27] + armWWbar[28];
   armWWbar[28]=1./12.*MMH;
   armWWbar[27]=armWWbar[27]*armWWbar[28];
   armWWbar[19]=armWWbar[21]*armWWbar[19];
   armWWbar[19]=armWWbar[25] - armWWbar[27] + armWWbar[19];
   armWWbar[21]= - armWWbar[22] - 1;
   armWWbar[21]=armWWbar[22]*armWWbar[21];
   armWWbar[25]=pow(armWWbar[5],2);
   armWWbar[21]=armWWbar[21] - armWWbar[25];
   armWWbar[19]=armWWbar[4]*armWWbar[19]*armWWbar[21];
   armWWbar[17]= - 5./4.*armWWbar[15] - 13./12.*armWWbar[16] + 
   armWWbar[17];
   armWWbar[17]=armWWbar[2]*armWWbar[17];
   armWWbar[17]=armWWbar[17] - armWWbar[18];
   armWWbar[17]=armWWbar[17]*armWWbar[25];
   armWWbar[17]=armWWbar[19] + armWWbar[24] + armWWbar[17];
   armWWbar[17]=armWWbar[4]*armWWbar[17];
   armWWbar[18]= - armWWbar[7]*armWWbar[1];
   armWWbar[19]=1./2.*armWWbar[15] + armWWbar[16] + MMZ;
   armWWbar[19]=armWWbar[9]*armWWbar[19];
   armWWbar[21]= - 1./3. + armWWbar[7];
   armWWbar[21]=armWWbar[11]*armWWbar[21];
   armWWbar[18]=armWWbar[26] + 4./3.*armWWbar[21] - 33./4.*armWWbar[13]
    + 3*armWWbar[19] - 59./18. + armWWbar[18];
   armWWbar[18]=armWWbar[2]*armWWbar[18];
   armWWbar[18]=armWWbar[18] + armWWbar[20];
   armWWbar[18]=armWWbar[18]*armWWbar[25];
   armWWbar[19]=armWWbar[23] + 17./12.;
   armWWbar[19]=armWWbar[13]*armWWbar[19];
   armWWbar[20]=MMZ + 3./2.*armWWbar[15];
   armWWbar[20]=armWWbar[9]*armWWbar[20];
   armWWbar[19]= - 1./6. + armWWbar[20] + armWWbar[19];
   armWWbar[19]=armWWbar[22]*armWWbar[19];
   armWWbar[20]= - armWWbar[9]*MMZ;
   armWWbar[20]=2*armWWbar[13] - 2 + armWWbar[20];
   armWWbar[19]=2*armWWbar[20] + armWWbar[19];
   armWWbar[19]=armWWbar[2]*armWWbar[19];

      mWWbarret = armWWbar[17] + armWWbar[18] + armWWbar[19];
      return mWWbarret;
}
