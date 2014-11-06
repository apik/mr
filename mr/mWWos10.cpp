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
WW<MS>::x10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWos[27], mWWosret;

    armWWos[1]=double(nL + nH);
    armWWos[2]=pow(s,-1);
    armWWos[3]=std::real(Tsil::B(0,0,mmW,mu2));
    armWWos[4]=double(nH);
    armWWos[5]=pow(mmZ,-1);
    armWWos[6]=pow(c,-1);
    armWWos[7]=Tsil::B(0,mmt,mmW,mu2);
    armWWos[8]=Tsil::A(mmt,mu2);
    armWWos[9]=pow(mmH,-1);
    armWWos[10]=double(nL);
    armWWos[11]=double(boson);
    armWWos[12]=Tsil::B(mmW,mmZ,mmW,mu2);
    armWWos[13]=Tsil::B(mmW,mmH,mmW,mu2);
    armWWos[14]=Tsil::A(mmW,mu2);
    armWWos[15]=Tsil::A(mmZ,mu2);
    armWWos[16]=Tsil::A(mmH,mu2);
   armWWos[17]=1./2.*armWWos[5];
   armWWos[18]=armWWos[13]*armWWos[5];
   armWWos[18]=armWWos[17] + armWWos[18];
   armWWos[19]=pow(armWWos[6],2);
   armWWos[20]=pow(armWWos[2],2);
   armWWos[21]=armWWos[19] + armWWos[20];
   armWWos[18]=armWWos[21]*armWWos[18];
   armWWos[22]=armWWos[19] + 1;
   armWWos[22]=armWWos[22]*armWWos[19];
   armWWos[22]=armWWos[22] + armWWos[20];
   armWWos[23]=1./4.*armWWos[22];
   armWWos[24]= - mmH*armWWos[13];
   armWWos[24]=armWWos[24] + armWWos[14] - armWWos[16];
   armWWos[25]=pow(armWWos[5],2);
   armWWos[23]=armWWos[25]*armWWos[23]*armWWos[24];
   armWWos[18]=armWWos[23] + armWWos[18];
   armWWos[18]=mmH*armWWos[18];
   armWWos[23]= - 3 - 1./3.*armWWos[19];
   armWWos[23]=armWWos[23]*armWWos[19];
   armWWos[23]=armWWos[23] + 5*armWWos[20];
   armWWos[23]=armWWos[15]*armWWos[23];
   armWWos[24]= - 35 + armWWos[19];
   armWWos[24]=armWWos[24]*armWWos[19];
   armWWos[24]=armWWos[24] + 13*armWWos[20];
   armWWos[24]=armWWos[14]*armWWos[24];
   armWWos[23]=1./12.*armWWos[24] + 1./4.*armWWos[23];
   armWWos[23]=armWWos[5]*armWWos[23];
   armWWos[24]= - armWWos[13] + 59./18.;
   armWWos[24]=armWWos[20]*armWWos[24];
   armWWos[17]= - armWWos[16]*armWWos[21]*armWWos[17];
   armWWos[26]= - 17 - armWWos[19];
   armWWos[26]=armWWos[26]*armWWos[19];
   armWWos[26]=33./4.*armWWos[20] - 4 + 1./12.*armWWos[26];
   armWWos[26]=armWWos[12]*armWWos[26];
   armWWos[17]=1./3.*armWWos[18] + armWWos[26] + armWWos[17] + 4 + 1./6.
   *armWWos[19] + armWWos[24] + armWWos[23];
   armWWos[17]=armWWos[11]*armWWos[17];
   armWWos[18]=armWWos[21]*armWWos[5];
   armWWos[22]=armWWos[22]*armWWos[25];
   armWWos[23]=1./2.*mmt;
   armWWos[24]=armWWos[23]*armWWos[22];
   armWWos[24]= - armWWos[18] + armWWos[24];
   armWWos[24]=armWWos[8]*armWWos[24];
   armWWos[22]=mmt*armWWos[22];
   armWWos[22]=armWWos[18] + armWWos[22];
   armWWos[22]=armWWos[22]*armWWos[23];
   armWWos[22]= - armWWos[20] + armWWos[22];
   armWWos[22]=armWWos[7]*armWWos[22];
   armWWos[18]=mmt*armWWos[18];
   armWWos[22]=armWWos[22] + armWWos[24] + 1./3.*armWWos[20] - 
   armWWos[18];
   armWWos[22]=armWWos[4]*armWWos[22];
   armWWos[23]= - 1./3.*armWWos[1] - armWWos[10];
   armWWos[24]=armWWos[3] - 1./3.;
   armWWos[23]=armWWos[23]*armWWos[24]*armWWos[20];
   armWWos[20]=3*armWWos[20];
   armWWos[19]= - armWWos[20] + 2 - armWWos[19];
   armWWos[19]=mmZ*armWWos[19];
   armWWos[21]=armWWos[15]*armWWos[21];
   armWWos[20]= - armWWos[14]*armWWos[20];
   armWWos[19]=armWWos[20] + armWWos[19] - 3./2.*armWWos[21];
   armWWos[19]=armWWos[11]*armWWos[19];
   armWWos[18]=armWWos[4]*armWWos[8]*armWWos[18];
   armWWos[18]=6*armWWos[18] + armWWos[19];
   armWWos[18]=armWWos[9]*armWWos[18];

      mWWosret = armWWos[17] + armWWos[18] + armWWos[22] + armWWos[23];
      return mWWosret;
}
