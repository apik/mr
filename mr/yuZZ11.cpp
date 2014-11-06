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
ZZ<OS>::y11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuZZ[32], yuZZret;

    aryuZZ[1]=double(nH);
    aryuZZ[2]=pow(CW,-1);
    aryuZZ[3]=pow(MMZ,-1);
    aryuZZ[4]=pow(SW,-1);
    aryuZZ[5]=Tsil::I2(0,0,MMt,mu2);
    aryuZZ[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuZZ[7]=Tsil::A(MMt,mu2);
    aryuZZ[8]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    aryuZZ[9]=pow(MMt,-1);
    aryuZZ[10]=Tsil::Aeps(MMt,mu2);
    aryuZZ[11]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryuZZ[12]=prottttt0->M(0);
    aryuZZ[13]=prot00000->M(0);
    aryuZZ[14]=prottttt0->Vzxyv(0);
    aryuZZ[15]=prottttt0->Suxv(0);
    aryuZZ[16]=double(nL);
    aryuZZ[17]=1/(4*MMt - MMZ);
   aryuZZ[18]=pow(aryuZZ[4],2);
   aryuZZ[19]=pow(aryuZZ[2],2);
   aryuZZ[20]=aryuZZ[18] + aryuZZ[19];
   aryuZZ[21]=aryuZZ[20]*aryuZZ[3];
   aryuZZ[22]=4*aryuZZ[7] - aryuZZ[15];
   aryuZZ[22]=aryuZZ[22]*aryuZZ[21];
   aryuZZ[22]=4*aryuZZ[22] - 197./18.*aryuZZ[19] - 128./9. - 29./2.*
   aryuZZ[18];
   aryuZZ[22]=aryuZZ[3]*aryuZZ[22];
   aryuZZ[23]=aryuZZ[20]*aryuZZ[7];
   aryuZZ[24]=2*aryuZZ[3];
   aryuZZ[25]=aryuZZ[24]*aryuZZ[23];
   aryuZZ[25]=aryuZZ[25] + 29./3.*aryuZZ[19] - 128./3. - aryuZZ[18];
   aryuZZ[25]=aryuZZ[25]*aryuZZ[24];
   aryuZZ[26]= - 7./9.*aryuZZ[19] + aryuZZ[18] + 64./9.;
   aryuZZ[27]=aryuZZ[26]*aryuZZ[3];
   aryuZZ[28]= - aryuZZ[6]*aryuZZ[27];
   aryuZZ[25]=aryuZZ[25] + aryuZZ[28];
   aryuZZ[25]=aryuZZ[6]*aryuZZ[25];
   aryuZZ[26]=aryuZZ[12]*aryuZZ[26];
   aryuZZ[26]=2*aryuZZ[26] - aryuZZ[21];
   aryuZZ[26]=aryuZZ[3]*aryuZZ[26];
   aryuZZ[28]=aryuZZ[20]*pow(aryuZZ[3],2);
   aryuZZ[29]=2 + aryuZZ[6];
   aryuZZ[29]=aryuZZ[6]*aryuZZ[29]*aryuZZ[28];
   aryuZZ[26]=aryuZZ[26] + aryuZZ[29];
   aryuZZ[26]=MMt*aryuZZ[26];
   aryuZZ[29]=MMt*aryuZZ[27];
   aryuZZ[30]=17./9.*aryuZZ[19] + aryuZZ[18] - 32./9.;
   aryuZZ[31]=aryuZZ[29] - aryuZZ[30];
   aryuZZ[31]=aryuZZ[14]*aryuZZ[31];
   aryuZZ[20]=aryuZZ[12]*aryuZZ[20];
   aryuZZ[20]=16./3.*aryuZZ[31] + 4./3.*aryuZZ[26] + 2./3.*aryuZZ[25]
    - 4*aryuZZ[20] + 1./3.*aryuZZ[22];
   aryuZZ[20]=MMt*aryuZZ[20];
   aryuZZ[22]=25./9.*aryuZZ[19] + aryuZZ[18] - 64./9.;
   aryuZZ[25]=aryuZZ[22]*aryuZZ[6];
   aryuZZ[26]=aryuZZ[7]*aryuZZ[27];
   aryuZZ[26]= - aryuZZ[25] - 8./3.*aryuZZ[26] + 143./9.*aryuZZ[19] - 
   320./9. + 7*aryuZZ[18];
   aryuZZ[26]=aryuZZ[6]*aryuZZ[26];
   aryuZZ[27]=41./9.*aryuZZ[19] - 128./9. + aryuZZ[18];
   aryuZZ[23]=aryuZZ[9]*aryuZZ[23];
   aryuZZ[23]=7./3.*aryuZZ[27] - 2*aryuZZ[23];
   aryuZZ[23]=aryuZZ[7]*aryuZZ[23];
   aryuZZ[27]= - aryuZZ[15]*aryuZZ[30];
   aryuZZ[23]=2*aryuZZ[27] + aryuZZ[23];
   aryuZZ[23]=aryuZZ[23]*aryuZZ[24];
   aryuZZ[27]=4./3.*aryuZZ[29] - aryuZZ[22];
   aryuZZ[27]=aryuZZ[27]*aryuZZ[8];
   aryuZZ[24]=aryuZZ[30]*aryuZZ[24];
   aryuZZ[28]=MMt*aryuZZ[28];
   aryuZZ[24]=aryuZZ[24] + aryuZZ[28];
   aryuZZ[24]=aryuZZ[24]*aryuZZ[10];
   aryuZZ[28]=937./18.*aryuZZ[19] - 860./9. + 77./2.*aryuZZ[18];
   aryuZZ[20]=8./3.*aryuZZ[24] + aryuZZ[27] + aryuZZ[26] + 1./3.*
   aryuZZ[28] + aryuZZ[23] + aryuZZ[20];
   aryuZZ[20]=aryuZZ[1]*aryuZZ[20];
   aryuZZ[23]= - aryuZZ[25] + 25./3.*aryuZZ[19] - 64./3. + 3*aryuZZ[18]
   ;
   aryuZZ[23]=aryuZZ[6]*aryuZZ[23];
   aryuZZ[24]= - aryuZZ[8]*aryuZZ[22];
   aryuZZ[23]=aryuZZ[24] + aryuZZ[23];
   aryuZZ[24]=MMZ*aryuZZ[1];
   aryuZZ[23]=aryuZZ[24]*aryuZZ[23];
   aryuZZ[25]=1 - aryuZZ[6];
   aryuZZ[25]=aryuZZ[7]*aryuZZ[25];
   aryuZZ[25]= - aryuZZ[10] + aryuZZ[25];
   aryuZZ[26]=4*aryuZZ[1];
   aryuZZ[22]=aryuZZ[26]*aryuZZ[22]*aryuZZ[25];
   aryuZZ[22]=aryuZZ[22] + aryuZZ[23];
   aryuZZ[22]=aryuZZ[17]*aryuZZ[22];
   aryuZZ[23]=aryuZZ[18] - 8./9. + 5./9.*aryuZZ[19];
   aryuZZ[23]=aryuZZ[23]*aryuZZ[1];
   aryuZZ[18]=11./9.*aryuZZ[19] + aryuZZ[18] - 20./9.;
   aryuZZ[19]=2*aryuZZ[16];
   aryuZZ[19]=aryuZZ[19]*aryuZZ[18];
   aryuZZ[19]=aryuZZ[23] + aryuZZ[19];
   aryuZZ[23]=aryuZZ[13]*MMZ;
   aryuZZ[23]=2*aryuZZ[11] + 4./3.*aryuZZ[23];
   aryuZZ[19]=aryuZZ[19]*aryuZZ[23];
   aryuZZ[18]=aryuZZ[16]*aryuZZ[18];
   aryuZZ[23]=aryuZZ[12]*aryuZZ[30]*aryuZZ[24];
   aryuZZ[21]=aryuZZ[5]*aryuZZ[26]*aryuZZ[21];

      yuZZret = 31./3.*aryuZZ[18] + aryuZZ[19] + aryuZZ[20] + 
      aryuZZ[21] + aryuZZ[22] + 4./3.*aryuZZ[23];
      return yuZZret;
}
