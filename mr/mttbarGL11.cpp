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
tt::mgl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[30], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=pow(SW,-1);
    armttbarGL[3]=pow(MMH,-1);
    armttbarGL[4]=pow(MMW,-1);
    armttbarGL[5]=pow(MMt,-1);
    armttbarGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    armttbarGL[7]=Tsil::I2(0,0,MMt,mu2);
    armttbarGL[8]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbarGL[9]=Tsil::A(MMH,mu2);
    armttbarGL[10]=Tsil::A(MMt,mu2);
    armttbarGL[11]=Tsil::Beps(MMH,MMt,MMt,mu2);
    armttbarGL[12]=std::real(Tsil::B(0,0,MMt,mu2));
    armttbarGL[13]=Tsil::Aeps(MMH,mu2);
    armttbarGL[14]=Tsil::Aeps(MMt,mu2);
    armttbarGL[15]=prot0ttHt->M(0);
    armttbarGL[16]=prottH0H->Vxzuv(0);
    armttbarGL[17]=prot0ttHt->Tuxv(0);
    armttbarGL[18]=protWt000->Tyzv(0);
    armttbarGL[19]=protHtt0t->M(0);
    armttbarGL[20]=prot00t00->M(0);
    armttbarGL[21]=prot000t0->M(0);
   armttbarGL[22]=2*MMH;
   armttbarGL[23]=armttbarGL[15] + armttbarGL[19];
   armttbarGL[24]= - 8./3.*armttbarGL[16] + armttbarGL[23];
   armttbarGL[24]=armttbarGL[24]*armttbarGL[22];
   armttbarGL[25]=armttbarGL[20] + armttbarGL[21];
   armttbarGL[25]=32*armttbarGL[3] - 8./3.*armttbarGL[23] - 1./3.*
   armttbarGL[25];
   armttbarGL[25]=MMt*armttbarGL[25];
   armttbarGL[26]=52*armttbarGL[8] - 17./2. - 8*armttbarGL[11];
   armttbarGL[27]=13./3.*armttbarGL[13] + 4*armttbarGL[10];
   armttbarGL[27]=armttbarGL[3]*armttbarGL[27];
   armttbarGL[28]=1 + 1./2.*armttbarGL[12];
   armttbarGL[28]=armttbarGL[12]*armttbarGL[28];
   armttbarGL[29]=pow(Pi,2);
   armttbarGL[24]=armttbarGL[25] + 13./3.*armttbarGL[17] - 5./12.*
   armttbarGL[29] + armttbarGL[24] + 5./6.*armttbarGL[28] + 1./3.*
   armttbarGL[26] + armttbarGL[27];
   armttbarGL[24]=MMt*armttbarGL[24];
   armttbarGL[25]=armttbarGL[5]*MMH;
   armttbarGL[25]=armttbarGL[25] - 2;
   armttbarGL[26]=pow(armttbarGL[10],2);
   armttbarGL[25]=armttbarGL[26]*armttbarGL[25];
   armttbarGL[27]=1 + armttbarGL[8];
   armttbarGL[27]=MMH*armttbarGL[27];
   armttbarGL[27]=armttbarGL[6] + armttbarGL[27] - armttbarGL[9];
   armttbarGL[28]= - 1 - 7./2.*armttbarGL[8];
   armttbarGL[28]=armttbarGL[10]*armttbarGL[28];
   armttbarGL[27]= - 5*armttbarGL[14] + armttbarGL[28] + 1./2.*
   armttbarGL[27];
   armttbarGL[28]=1./3.*MMH;
   armttbarGL[27]=armttbarGL[27]*armttbarGL[28];
   armttbarGL[29]=armttbarGL[17]*pow(MMH,2);
   armttbarGL[25]=1./6.*armttbarGL[29] + armttbarGL[27] + 
   armttbarGL[25];
   armttbarGL[25]=armttbarGL[5]*armttbarGL[25];
   armttbarGL[23]=4*armttbarGL[16] - armttbarGL[23];
   armttbarGL[23]=MMH*armttbarGL[23];
   armttbarGL[23]=armttbarGL[23] - 17*armttbarGL[8] - 13 + 4*
   armttbarGL[11];
   armttbarGL[23]=armttbarGL[23]*armttbarGL[28];
   armttbarGL[27]=4./3.*armttbarGL[12] + 23./2. + 14./3.*armttbarGL[8];
   armttbarGL[27]=armttbarGL[10]*armttbarGL[27];
   armttbarGL[26]=armttbarGL[3]*armttbarGL[26];
   armttbarGL[29]= - 4./3.*armttbarGL[9] - 13./3.*armttbarGL[10];
   armttbarGL[29]=armttbarGL[3]*armttbarGL[29];
   armttbarGL[29]=7./2. - 2./3.*armttbarGL[8] + armttbarGL[29];
   armttbarGL[29]=armttbarGL[9]*armttbarGL[29];
   armttbarGL[22]= - armttbarGL[17]*armttbarGL[22];
   armttbarGL[28]= - armttbarGL[28] + MMt;
   armttbarGL[28]=armttbarGL[18]*armttbarGL[28];
   armttbarGL[22]=armttbarGL[25] + 4*armttbarGL[28] + armttbarGL[24] + 
   armttbarGL[22] + armttbarGL[23] + 65./6.*armttbarGL[14] - 13./6.*
   armttbarGL[7] + armttbarGL[29] - 36*armttbarGL[26] - 11./3.*
   armttbarGL[6] + 3*armttbarGL[13] + armttbarGL[27];

      mttbarGLret = armttbarGL[22]*armttbarGL[4]*pow(armttbarGL[2],2)*
      armttbarGL[1];
      return mttbarGLret;
}
