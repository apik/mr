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
tt::mygl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[31], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=pow(SW,-1);
    aryuttGL[3]=pow(MMW,-1);
    aryuttGL[4]=pow(MMt,-1);
    aryuttGL[5]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuttGL[6]=Tsil::I2(0,0,MMt,mu2);
    aryuttGL[7]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[8]=Tsil::A(MMH,mu2);
    aryuttGL[9]=Tsil::A(MMt,mu2);
    aryuttGL[10]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryuttGL[11]=pow(MMH,-1);
    aryuttGL[12]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuttGL[13]=Tsil::Aeps(MMH,mu2);
    aryuttGL[14]=Tsil::Aeps(MMt,mu2);
    aryuttGL[15]=prot0ttHt->M(0);
    aryuttGL[16]=prottH0H->Vxzuv(0);
    aryuttGL[17]=prot0ttHt->Tuxv(0);
    aryuttGL[18]=protWt000->Tyzv(0);
    aryuttGL[19]=protHtt0t->M(0);
    aryuttGL[20]=prot00t00->M(0);
    aryuttGL[21]=prot000t0->M(0);
   aryuttGL[22]=4*aryuttGL[18];
   aryuttGL[23]=aryuttGL[14]*aryuttGL[4];
   aryuttGL[24]=aryuttGL[15] + aryuttGL[19];
   aryuttGL[25]=4*aryuttGL[16] - aryuttGL[24] + 1./2.*aryuttGL[4];
   aryuttGL[25]=MMH*aryuttGL[25];
   aryuttGL[26]=MMH*aryuttGL[4];
   aryuttGL[27]=1./2.*aryuttGL[26];
   aryuttGL[28]= - 17 + aryuttGL[27];
   aryuttGL[28]=aryuttGL[7]*aryuttGL[28];
   aryuttGL[23]=aryuttGL[28] + aryuttGL[25] - 5*aryuttGL[23] + 4*
   aryuttGL[10] - 53./4. - aryuttGL[22];
   aryuttGL[25]=1./3.*MMH;
   aryuttGL[23]=aryuttGL[25]*aryuttGL[23];
   aryuttGL[25]=2 + 1./3.*aryuttGL[12];
   aryuttGL[28]=2 - aryuttGL[27];
   aryuttGL[28]=aryuttGL[7]*aryuttGL[28];
   aryuttGL[29]=MMH*pow(aryuttGL[4],2);
   aryuttGL[29]= - 8*aryuttGL[4] + aryuttGL[29];
   aryuttGL[29]=aryuttGL[9]*aryuttGL[29];
   aryuttGL[25]=aryuttGL[29] + 7./3.*aryuttGL[28] + 4*aryuttGL[25] - 1./
   12.*aryuttGL[26];
   aryuttGL[25]=aryuttGL[9]*aryuttGL[25];
   aryuttGL[26]=1./6.*aryuttGL[26];
   aryuttGL[28]=13./3.*aryuttGL[11];
   aryuttGL[29]= - 3./2.*aryuttGL[4] - aryuttGL[28];
   aryuttGL[29]=aryuttGL[9]*aryuttGL[29];
   aryuttGL[30]=aryuttGL[8]*aryuttGL[11];
   aryuttGL[29]= - 4./3.*aryuttGL[30] + aryuttGL[29] - 2./3.*
   aryuttGL[7] + 4 - aryuttGL[26];
   aryuttGL[29]=aryuttGL[8]*aryuttGL[29];
   aryuttGL[27]= - 11 + aryuttGL[27];
   aryuttGL[27]=aryuttGL[5]*aryuttGL[27];
   aryuttGL[23]=aryuttGL[29] + 1./3.*aryuttGL[27] + aryuttGL[25] - 13./
   6.*aryuttGL[6] + 65./6.*aryuttGL[14] + 3*aryuttGL[13] + aryuttGL[23]
   ;
   aryuttGL[23]=aryuttGL[3]*aryuttGL[23];
   aryuttGL[25]=1 + 1./2.*aryuttGL[12];
   aryuttGL[25]=aryuttGL[12]*aryuttGL[25];
   aryuttGL[27]=aryuttGL[13]*aryuttGL[28];
   aryuttGL[28]= - 8./3.*aryuttGL[16] + aryuttGL[24];
   aryuttGL[28]=MMH*aryuttGL[28];
   aryuttGL[29]=pow(Pi,2);
   aryuttGL[22]=52./3.*aryuttGL[7] - 3./4.*aryuttGL[29] + 2*
   aryuttGL[28] + aryuttGL[27] - 8./3.*aryuttGL[10] + 5./6.*
   aryuttGL[25] + 59./12. + aryuttGL[22];
   aryuttGL[22]=aryuttGL[3]*aryuttGL[22];
   aryuttGL[25]=MMt*aryuttGL[3];
   aryuttGL[24]= - aryuttGL[20] - aryuttGL[21] - 8*aryuttGL[24];
   aryuttGL[24]=aryuttGL[24]*aryuttGL[25];
   aryuttGL[22]=aryuttGL[22] + 1./3.*aryuttGL[24];
   aryuttGL[22]=MMt*aryuttGL[22];
   aryuttGL[24]= - 2 + aryuttGL[26];
   aryuttGL[24]=aryuttGL[3]*MMH*aryuttGL[24];
   aryuttGL[24]=aryuttGL[24] + 13./3.*aryuttGL[25];
   aryuttGL[24]=aryuttGL[17]*aryuttGL[24];
   aryuttGL[22]=aryuttGL[24] + aryuttGL[23] + aryuttGL[22];

      yuttGLret = aryuttGL[22]*pow(aryuttGL[2],2)*aryuttGL[1];
      return yuttGLret;
}
