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
bb::mgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbarGL[28], mbbbarGLret;

    armbbbarGL[1]=double(boson);
    armbbbarGL[2]=pow(SW,-1);
    armbbbarGL[3]=pow(MMH,-1);
    armbbbarGL[4]=pow(MMW,-1);
    armbbbarGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    armbbbarGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    armbbbarGL[7]=pow(MMt,-1);
    armbbbarGL[8]=Tsil::I2(0,MMH,MMt,mu2);
    armbbbarGL[9]=Tsil::I2(0,0,MMH,mu2);
    armbbbarGL[10]=Tsil::I2(0,0,MMt,mu2);
    armbbbarGL[11]=Tsil::B(MMH,MMH,MMH,mu2);
    armbbbarGL[12]=Tsil::A(MMH,mu2);
    armbbbarGL[13]=Tsil::A(MMt,mu2);
    armbbbarGL[14]=Tsil::B(MMH,MMt,MMt,mu2);
    armbbbarGL[15]=Tsil::B(MMt,MMt,MMH,mu2);
    armbbbarGL[16]=std::real(Tsil::B(0,0,MMH,mu2));
    armbbbarGL[17]=std::real(Tsil::B(0,0,MMt,mu2));
    armbbbarGL[18]=Tsil::Aeps(MMH,mu2);
    armbbbarGL[19]=Tsil::Aeps(MMt,mu2);
   armbbbarGL[20]=15./2.*armbbbarGL[18] - 7*armbbbarGL[6] - 1./2.*
   armbbbarGL[8];
   armbbbarGL[21]=9*armbbbarGL[15];
   armbbbarGL[22]=49./8. + armbbbarGL[21];
   armbbbarGL[22]=armbbbarGL[12]*armbbbarGL[22];
   armbbbarGL[23]= - 85./8. + armbbbarGL[21];
   armbbbarGL[23]=MMH*armbbbarGL[23];
   armbbbarGL[20]=armbbbarGL[23] - 31./4.*armbbbarGL[10] + 5./2.*
   armbbbarGL[20] + armbbbarGL[22];
   armbbbarGL[22]=3*armbbbarGL[3];
   armbbbarGL[23]=3*armbbbarGL[15];
   armbbbarGL[24]= - 5./2. - armbbbarGL[23];
   armbbbarGL[22]=armbbbarGL[12]*armbbbarGL[24]*armbbbarGL[22];
   armbbbarGL[24]=27*armbbbarGL[6] - armbbbarGL[8];
   armbbbarGL[24]=13./2.*armbbbarGL[10] + 1./2.*armbbbarGL[24] - 13*
   armbbbarGL[18];
   armbbbarGL[24]=armbbbarGL[3]*armbbbarGL[24];
   armbbbarGL[21]=armbbbarGL[22] + 7./16.*armbbbarGL[17] + 463./128. - 
   armbbbarGL[21] + armbbbarGL[24];
   armbbbarGL[22]=MMt*armbbbarGL[3];
   armbbbarGL[24]=3*armbbbarGL[22];
   armbbbarGL[25]=3 - 1./2.*armbbbarGL[17];
   armbbbarGL[25]=armbbbarGL[25]*armbbbarGL[24];
   armbbbarGL[21]=1./2.*armbbbarGL[21] + armbbbarGL[25];
   armbbbarGL[21]=MMt*armbbbarGL[21];
   armbbbarGL[20]=1./8.*armbbbarGL[20] + armbbbarGL[21];
   armbbbarGL[20]=MMt*armbbbarGL[20];
   armbbbarGL[21]=19*armbbbarGL[8] - 17*armbbbarGL[9] - 13*
   armbbbarGL[6];
   armbbbarGL[21]= - 15./8.*armbbbarGL[5] + 31./8.*armbbbarGL[12] + 1./
   8.*armbbbarGL[21] + 7*armbbbarGL[18];
   armbbbarGL[25]= - 3./16.*armbbbarGL[8] + 1./16.*armbbbarGL[6] + 1./8.
   *armbbbarGL[9];
   armbbbarGL[25]=armbbbarGL[7]*armbbbarGL[25];
   armbbbarGL[25]= - 1 + armbbbarGL[25];
   armbbbarGL[25]=MMH*armbbbarGL[25];
   armbbbarGL[21]=1./4.*armbbbarGL[21] + armbbbarGL[25];
   armbbbarGL[21]=MMH*armbbbarGL[21];
   armbbbarGL[25]=pow(armbbbarGL[12],2);
   armbbbarGL[20]=armbbbarGL[20] - 9./64.*armbbbarGL[25] + 
   armbbbarGL[21];
   armbbbarGL[21]=61./4. + 3*armbbbarGL[17];
   armbbbarGL[25]=armbbbarGL[12]*armbbbarGL[3];
   armbbbarGL[21]=1./4.*armbbbarGL[21] + 11*armbbbarGL[25];
   armbbbarGL[25]=1./2.*armbbbarGL[3];
   armbbbarGL[23]= - armbbbarGL[17] - 13./8. + armbbbarGL[23];
   armbbbarGL[23]=armbbbarGL[23]*armbbbarGL[25];
   armbbbarGL[26]=MMt*pow(armbbbarGL[3],2);
   armbbbarGL[27]=armbbbarGL[15]*armbbbarGL[26];
   armbbbarGL[23]=armbbbarGL[23] - 6*armbbbarGL[27];
   armbbbarGL[27]=3*MMt;
   armbbbarGL[23]=armbbbarGL[23]*armbbbarGL[27];
   armbbbarGL[21]=1./8.*armbbbarGL[21] + armbbbarGL[23];
   armbbbarGL[21]=MMt*armbbbarGL[21];
   armbbbarGL[23]=armbbbarGL[25] - armbbbarGL[26];
   armbbbarGL[23]=MMt*armbbbarGL[23];
   armbbbarGL[25]=MMH*armbbbarGL[7];
   armbbbarGL[26]=63./4. - armbbbarGL[25];
   armbbbarGL[23]=1./16.*armbbbarGL[26] + 9*armbbbarGL[23];
   armbbbarGL[23]=armbbbarGL[13]*armbbbarGL[23];
   armbbbarGL[26]= - armbbbarGL[12]*armbbbarGL[7];
   armbbbarGL[26]=1./2. + armbbbarGL[26];
   armbbbarGL[26]=MMH*armbbbarGL[26];
   armbbbarGL[26]= - 9*armbbbarGL[12] + 1./2.*armbbbarGL[26];
   armbbbarGL[21]=1./2.*armbbbarGL[23] + 1./16.*armbbbarGL[26] + 
   armbbbarGL[21];
   armbbbarGL[21]=armbbbarGL[13]*armbbbarGL[21];
   armbbbarGL[23]=MMH + armbbbarGL[12];
   armbbbarGL[26]=1./8.*MMH;
   armbbbarGL[23]=armbbbarGL[23]*armbbbarGL[26];
   armbbbarGL[27]=armbbbarGL[13]*MMt;
   armbbbarGL[23]=armbbbarGL[23] + armbbbarGL[27];
   armbbbarGL[27]=27./8.*armbbbarGL[11] + 9./8.*armbbbarGL[16];
   armbbbarGL[23]=armbbbarGL[23]*armbbbarGL[27];
   armbbbarGL[25]=7./2. + armbbbarGL[25];
   armbbbarGL[25]=armbbbarGL[25]*armbbbarGL[26];
   armbbbarGL[26]=1 - armbbbarGL[24];
   armbbbarGL[26]=MMt*armbbbarGL[26];
   armbbbarGL[25]=armbbbarGL[25] + 11*armbbbarGL[26];
   armbbbarGL[25]=armbbbarGL[19]*armbbbarGL[25];
   armbbbarGL[24]=19./16. - armbbbarGL[24];
   armbbbarGL[24]=MMt*armbbbarGL[24];
   armbbbarGL[24]= - 7./64.*MMH + armbbbarGL[24];
   armbbbarGL[24]=MMt*armbbbarGL[24];
   armbbbarGL[22]=5./8. - 2*armbbbarGL[22];
   armbbbarGL[22]=MMt*armbbbarGL[22];
   armbbbarGL[22]= - 1./32.*MMH + armbbbarGL[22];
   armbbbarGL[22]=armbbbarGL[13]*armbbbarGL[22];
   armbbbarGL[22]=armbbbarGL[24] + 3*armbbbarGL[22];
   armbbbarGL[22]=armbbbarGL[14]*armbbbarGL[22];
   armbbbarGL[20]=armbbbarGL[22] + armbbbarGL[23] + 1./4.*
   armbbbarGL[25] + 1./2.*armbbbarGL[20] + armbbbarGL[21];

      mbbbarGLret = armbbbarGL[20]*pow(armbbbarGL[4],2)*pow(
      armbbbarGL[2],4)*armbbbarGL[1];
      return mbbbarGLret;
}
