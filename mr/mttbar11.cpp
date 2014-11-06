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
tt::x11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbar[75], mttbarret;

    armttbar[1]=double(nH);
    armttbar[2]=pow(CW,-1);
    armttbar[3]=pow(MMH,-1);
    armttbar[4]=pow(MMZ,-1);
    armttbar[5]=pow(SW,-1);
    armttbar[6]=Tsil::A(MMt,mu2);
    armttbar[7]=double(boson);
    armttbar[8]=pow(MMt,-1);
    armttbar[9]=Tsil::I2(MMH,MMt,MMt,mu2);
    armttbar[10]=Tsil::I2(MMZ,MMt,MMt,mu2);
    armttbar[11]=Tsil::I2(0,MMW,MMt,mu2);
    armttbar[12]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbar[13]=Tsil::A(MMH,mu2);
    armttbar[14]=Tsil::B(MMZ,MMt,MMt,mu2);
    armttbar[15]=Tsil::A(MMZ,mu2);
    armttbar[16]=Tsil::Beps(MMH,MMt,MMt,mu2);
    armttbar[17]=Tsil::Beps(MMZ,MMt,MMt,mu2);
    armttbar[18]=Tsil::A(MMW,mu2);
    armttbar[19]=std::real(Tsil::B(0,MMW,MMt,mu2));
    armttbar[20]=Tsil::Aeps(MMH,mu2);
    armttbar[21]=Tsil::Aeps(MMZ,mu2);
    armttbar[22]=Tsil::Aeps(MMW,mu2);
    armttbar[23]=Tsil::Aeps(MMt,mu2);
    armttbar[24]=std::real(Tsil::Beps(0,MMW,MMt,mu2));
    armttbar[25]=protWt000->M(0);
    armttbar[26]=prot0ttHt->M(0);
    armttbar[27]=prot0ttZt->M(0);
    armttbar[28]=prot0tt0t->M(0);
    armttbar[29]=prottH0H->Vxzuv(0);
    armttbar[30]=prottZ0Z->Vxzuv(0);
    armttbar[31]=protWt000->Uzxyv(0);
    armttbar[32]=prot0ttHt->Tuxv(0);
    armttbar[33]=prot0ttZt->Tuxv(0);
    armttbar[34]=protWt000->Txuv(0);
    armttbar[35]=protWt000->Tyzv(0);
    armttbar[36]=1/(4*MMt - MMZ);
   armttbar[37]=pow(armttbar[2],2);
   armttbar[38]=17./9.*armttbar[37];
   armttbar[39]=pow(armttbar[5],2);
   armttbar[40]=armttbar[38] + armttbar[39] - 32./9.;
   armttbar[41]= - armttbar[27] + 2*armttbar[30];
   armttbar[41]=armttbar[40]*armttbar[41];
   armttbar[42]=armttbar[39] - 1;
   armttbar[43]=pow(CW,2);
   armttbar[44]=armttbar[42] - armttbar[43];
   armttbar[45]=armttbar[25]*armttbar[44];
   armttbar[41]= - 2*armttbar[45] + armttbar[41];
   armttbar[41]=MMZ*armttbar[41];
   armttbar[45]=pow(Pi,2);
   armttbar[46]=5./2.*armttbar[45];
   armttbar[47]= - armttbar[46] + 2*armttbar[24];
   armttbar[48]= - 1 - armttbar[47];
   armttbar[48]=armttbar[48]*armttbar[39];
   armttbar[49]=95./18.*armttbar[37] - 64./9. + 7./2.*armttbar[39];
   armttbar[49]=1./3.*armttbar[49];
   armttbar[50]= - armttbar[14]*armttbar[49];
   armttbar[51]=3./2.*armttbar[39];
   armttbar[52]=59./18.*armttbar[37] - 64./9. + armttbar[51];
   armttbar[52]=armttbar[33]*armttbar[52];
   armttbar[53]=23 - 17*armttbar[37];
   armttbar[53]=armttbar[35]*armttbar[53];
   armttbar[54]=armttbar[34]*armttbar[42];
   armttbar[41]=armttbar[52] + 2./3.*armttbar[41] + armttbar[50] - 19./
   6.*armttbar[54] + 4./27.*armttbar[53] - 163./54.*armttbar[37] + 
   armttbar[48] + 211./54. + armttbar[47];
   armttbar[41]=MMZ*armttbar[41];
   armttbar[45]= - armttbar[34] + 5*armttbar[45];
   armttbar[47]= - armttbar[42]*armttbar[45];
   armttbar[48]=armttbar[40]*armttbar[14];
   armttbar[45]= - 1 + armttbar[45];
   armttbar[43]=armttbar[45]*armttbar[43];
   armttbar[43]=armttbar[43] - armttbar[48] - armttbar[38] + 23./9. + 
   armttbar[47];
   armttbar[43]=MMZ*armttbar[43];
   armttbar[45]=5*armttbar[39];
   armttbar[38]=armttbar[38] + 22./9. - armttbar[45];
   armttbar[47]=1./2.*armttbar[39];
   armttbar[50]=17./18.*armttbar[37] + armttbar[47] - 16./9.;
   armttbar[52]=armttbar[14]*armttbar[50];
   armttbar[38]= - 7*armttbar[52] + 2*armttbar[38];
   armttbar[38]=armttbar[6]*armttbar[38];
   armttbar[52]= - armttbar[18] + armttbar[11];
   armttbar[52]=armttbar[42]*armttbar[52];
   armttbar[53]= - 34./9.*armttbar[37] + 37./9. + armttbar[39];
   armttbar[53]=armttbar[23]*armttbar[53];
   armttbar[54]= - armttbar[10]*armttbar[40];
   armttbar[38]=armttbar[43] + armttbar[54] + armttbar[53] + 
   armttbar[38] + armttbar[52];
   armttbar[38]=MMZ*armttbar[38];
   armttbar[43]=armttbar[8]*MMZ;
   armttbar[52]=armttbar[40]*armttbar[43];
   armttbar[53]=armttbar[52] + armttbar[39] + 41./9.*armttbar[37];
   armttbar[54]=pow(armttbar[6],2);
   armttbar[53]=armttbar[54]*armttbar[53];
   armttbar[55]=armttbar[37] + armttbar[39];
   armttbar[56]=MMH*armttbar[55];
   armttbar[54]=armttbar[54]*armttbar[4];
   armttbar[57]=armttbar[56]*armttbar[54];
   armttbar[58]=1./3.*armttbar[40];
   armttbar[59]= - armttbar[33]*pow(MMZ,2)*armttbar[58];
   armttbar[60]=armttbar[39]*armttbar[18];
   armttbar[61]=armttbar[6]*armttbar[60];
   armttbar[38]=armttbar[59] + 1./3.*armttbar[38] - 5./3.*armttbar[61]
    + armttbar[57] + armttbar[53];
   armttbar[38]=armttbar[8]*armttbar[38];
   armttbar[53]=19./3.*armttbar[37] - 64./3. + armttbar[39];
   armttbar[53]=armttbar[6]*armttbar[53];
   armttbar[53]= - 1./3.*armttbar[56] + armttbar[53];
   armttbar[53]=armttbar[6]*armttbar[53];
   armttbar[57]=pow(MMH,2);
   armttbar[59]=armttbar[57]*armttbar[55];
   armttbar[61]=1./6.*armttbar[59];
   armttbar[62]=armttbar[23]*armttbar[56];
   armttbar[53]= - 5./3.*armttbar[62] + armttbar[61] + armttbar[53];
   armttbar[53]=armttbar[4]*armttbar[53];
   armttbar[62]= - 7./18.*armttbar[37] + 32./9. + armttbar[47];
   armttbar[62]=armttbar[14]*armttbar[62];
   armttbar[51]=7./3.*armttbar[62] + 119./54.*armttbar[37] + 320./27.
    + armttbar[51];
   armttbar[51]=armttbar[6]*armttbar[51];
   armttbar[62]=armttbar[6]*armttbar[56];
   armttbar[59]= - 7*armttbar[62] + armttbar[59];
   armttbar[62]=armttbar[12]*armttbar[4];
   armttbar[59]=armttbar[59]*armttbar[62];
   armttbar[63]=armttbar[55]*armttbar[4];
   armttbar[64]=armttbar[63]*armttbar[18];
   armttbar[65]=armttbar[4]*armttbar[6];
   armttbar[66]=armttbar[55]*armttbar[65];
   armttbar[67]=11./3.*armttbar[39] - 3*armttbar[66];
   armttbar[67]=1./2.*armttbar[67] + 5./3.*armttbar[64];
   armttbar[67]=armttbar[18]*armttbar[67];
   armttbar[68]= - 23./3.*armttbar[37] - 32./3. + 23./2.*armttbar[39];
   armttbar[68]=armttbar[23]*armttbar[68];
   armttbar[49]= - armttbar[10]*armttbar[49];
   armttbar[69]=armttbar[11]*armttbar[39];
   armttbar[38]=armttbar[38] + armttbar[49] + armttbar[67] + 1./6.*
   armttbar[59] + armttbar[53] + 1./3.*armttbar[68] - 13./6.*
   armttbar[69] + armttbar[41] + armttbar[51];
   armttbar[38]=armttbar[8]*armttbar[38];
   armttbar[41]=armttbar[29]*MMH;
   armttbar[41]=16*armttbar[41] + 8*armttbar[16] + armttbar[46] + 29./4.
   ;
   armttbar[46]=4*MMH;
   armttbar[49]=armttbar[46]*armttbar[26];
   armttbar[51]=2*armttbar[14];
   armttbar[41]=armttbar[49] - 1./3.*armttbar[41] - armttbar[51] + 5./6.
   *armttbar[34];
   armttbar[41]=armttbar[55]*armttbar[41];
   armttbar[49]=3*armttbar[55];
   armttbar[53]=armttbar[35]*armttbar[49];
   armttbar[41]=armttbar[53] + armttbar[41];
   armttbar[41]=armttbar[4]*armttbar[41];
   armttbar[53]=armttbar[37] + 1;
   armttbar[53]=armttbar[53]*armttbar[37];
   armttbar[53]=armttbar[53] + armttbar[39];
   armttbar[59]=pow(armttbar[4],2);
   armttbar[53]=armttbar[53]*armttbar[59];
   armttbar[67]=armttbar[18]*armttbar[53];
   armttbar[68]= - armttbar[33] + 52./3.*armttbar[12];
   armttbar[68]=armttbar[63]*armttbar[68];
   armttbar[69]= - 7./9.*armttbar[37] + armttbar[39] + 64./9.;
   armttbar[70]= - armttbar[27]*armttbar[69];
   armttbar[70]=64./9.*armttbar[28] + armttbar[70];
   armttbar[71]= - 8*armttbar[26] - armttbar[25];
   armttbar[71]=MMt*armttbar[63]*armttbar[71];
   armttbar[41]=2./3.*armttbar[71] - 5./6.*armttbar[67] + 4./3.*
   armttbar[70] + armttbar[41] + armttbar[68];
   armttbar[41]=MMt*armttbar[41];
   armttbar[50]=armttbar[6]*armttbar[50];
   armttbar[58]=MMZ*armttbar[58];
   armttbar[50]= - 5*armttbar[50] + armttbar[58];
   armttbar[50]=armttbar[8]*armttbar[50];
   armttbar[58]= - 35./9.*armttbar[37] + 104./9. - armttbar[39];
   armttbar[58]=armttbar[58]*armttbar[65];
   armttbar[48]=armttbar[50] + 4./3.*armttbar[58] - 2./3.*armttbar[48]
    + 77./18.*armttbar[37] - 64./9. + 5./2.*armttbar[39];
   armttbar[48]=armttbar[8]*armttbar[48];
   armttbar[50]= - 1./2.*armttbar[55] + armttbar[66];
   armttbar[50]=armttbar[4]*armttbar[50];
   armttbar[58]=armttbar[8]*armttbar[4];
   armttbar[68]=armttbar[15]*armttbar[40]*armttbar[58];
   armttbar[48]= - 4./3.*armttbar[68] + armttbar[50] + armttbar[48];
   armttbar[48]=armttbar[15]*armttbar[48];
   armttbar[50]=armttbar[15]*armttbar[4];
   armttbar[50]=armttbar[65] + armttbar[50];
   armttbar[50]= - 5 + armttbar[51] + 4*armttbar[50];
   armttbar[50]=armttbar[15]*armttbar[50];
   armttbar[68]=2*MMZ;
   armttbar[70]= - armttbar[17]*armttbar[68];
   armttbar[50]=armttbar[50] - 7*armttbar[21] + armttbar[70];
   armttbar[70]=25./9.*armttbar[37] + armttbar[39] - 64./9.;
   armttbar[50]=armttbar[70]*armttbar[50];
   armttbar[71]=3*armttbar[39];
   armttbar[72]=25./3.*armttbar[37] + armttbar[71] - 64./3.;
   armttbar[73]= - 2*armttbar[54] - armttbar[6];
   armttbar[73]=armttbar[72]*armttbar[73];
   armttbar[74]=armttbar[23]*armttbar[70];
   armttbar[73]=armttbar[74] + armttbar[73];
   armttbar[47]=25./18.*armttbar[37] - 32./9. + armttbar[47];
   armttbar[51]=armttbar[70]*armttbar[51];
   armttbar[47]=13*armttbar[47] + armttbar[51];
   armttbar[47]=MMZ*armttbar[47];
   armttbar[51]=armttbar[10]*armttbar[72];
   armttbar[47]=armttbar[47] + 2*armttbar[73] + armttbar[51] + 
   armttbar[50];
   armttbar[47]=armttbar[36]*armttbar[47];
   armttbar[44]=MMZ*armttbar[44];
   armttbar[50]=2*armttbar[42];
   armttbar[51]= - armttbar[6]*armttbar[50];
   armttbar[44]=armttbar[51] + armttbar[44];
   armttbar[44]=armttbar[44]*armttbar[43];
   armttbar[51]=armttbar[6]*armttbar[39];
   armttbar[70]=MMZ*armttbar[42];
   armttbar[44]=armttbar[44] - 47./2.*armttbar[70] + 4*armttbar[51] + 5
   *armttbar[60];
   armttbar[44]=armttbar[8]*armttbar[44];
   armttbar[44]=armttbar[44] - armttbar[39] + 4*armttbar[66];
   armttbar[51]=armttbar[63] + armttbar[67];
   armttbar[51]=MMt*armttbar[51];
   armttbar[70]=armttbar[63]*MMt;
   armttbar[72]=armttbar[39] + 1./2.*armttbar[70];
   armttbar[72]=armttbar[19]*armttbar[72];
   armttbar[44]=5./3.*armttbar[72] + 5./6.*armttbar[51] + 5./2.*
   armttbar[64] + 1./3.*armttbar[44];
   armttbar[44]=armttbar[19]*armttbar[44];
   armttbar[51]=armttbar[6]*armttbar[55];
   armttbar[51]= - 17*armttbar[56] + 14*armttbar[51];
   armttbar[51]=armttbar[51]*armttbar[62];
   armttbar[62]=armttbar[58]*armttbar[56];
   armttbar[64]= - 11*armttbar[63] + 1./2.*armttbar[62];
   armttbar[64]=armttbar[9]*armttbar[64];
   armttbar[51]=armttbar[64] + armttbar[51] - 47./24.*armttbar[37] - 
   704./9. - 119./8.*armttbar[39];
   armttbar[64]= - 2*armttbar[26] + 4*armttbar[29];
   armttbar[57]=armttbar[64]*armttbar[57];
   armttbar[46]=armttbar[46]*armttbar[16];
   armttbar[46]= - armttbar[57] - armttbar[46] + 13./2.*armttbar[11] + 
   13*MMH;
   armttbar[46]= - armttbar[46]*armttbar[55];
   armttbar[57]=armttbar[35]*armttbar[56];
   armttbar[46]= - 4*armttbar[57] + armttbar[46];
   armttbar[57]=59./6.*armttbar[23] + 23./2.*armttbar[6];
   armttbar[57]=armttbar[55]*armttbar[57];
   armttbar[64]=armttbar[20]*armttbar[49];
   armttbar[46]=armttbar[64] + 1./3.*armttbar[46] + armttbar[57];
   armttbar[46]=armttbar[4]*armttbar[46];
   armttbar[57]= - 1./3.*armttbar[37] + 2./3. - armttbar[39];
   armttbar[57]=MMZ*armttbar[57];
   armttbar[64]=armttbar[37] - 2 + armttbar[71];
   armttbar[64]=MMZ*armttbar[64];
   armttbar[71]=armttbar[18]*armttbar[71];
   armttbar[64]=armttbar[71] + armttbar[64];
   armttbar[71]=armttbar[8]*armttbar[6];
   armttbar[64]=armttbar[64]*armttbar[71];
   armttbar[57]=armttbar[64] - armttbar[60] + armttbar[57];
   armttbar[49]=armttbar[49]*armttbar[71];
   armttbar[49]=armttbar[49] - armttbar[55];
   armttbar[49]=armttbar[15]*armttbar[49];
   armttbar[60]=13./3.*armttbar[70];
   armttbar[64]=armttbar[20]*armttbar[60];
   armttbar[49]=armttbar[49] + armttbar[64] + 2*armttbar[57];
   armttbar[49]=armttbar[3]*armttbar[49];
   armttbar[57]=65./9.*armttbar[37] + armttbar[39] - 224./9.;
   armttbar[57]=1./3.*armttbar[57];
   armttbar[64]=armttbar[4]*armttbar[57];
   armttbar[59]= - armttbar[55]*armttbar[59]*MMt;
   armttbar[40]=armttbar[8]*armttbar[40];
   armttbar[40]=4*armttbar[40] + armttbar[64] + armttbar[59];
   armttbar[40]=armttbar[21]*armttbar[40];
   armttbar[59]=armttbar[13]*armttbar[3];
   armttbar[59]= - 4./3.*armttbar[59] + 7./2. - 2./3.*armttbar[12];
   armttbar[59]=armttbar[63]*armttbar[59];
   armttbar[64]=armttbar[3]*armttbar[66];
   armttbar[59]= - 13./3.*armttbar[64] - 1./6.*armttbar[62] + 
   armttbar[59];
   armttbar[59]=armttbar[13]*armttbar[59];
   armttbar[56]=armttbar[4]*armttbar[56];
   armttbar[58]=armttbar[58]*armttbar[61];
   armttbar[56]=armttbar[58] - 2*armttbar[56] + armttbar[60];
   armttbar[56]=armttbar[32]*armttbar[56];
   armttbar[52]=2*armttbar[52] - armttbar[69];
   armttbar[52]=armttbar[17]*armttbar[52];
   armttbar[45]=1./9.*armttbar[37] + 176./9. + armttbar[45];
   armttbar[45]=armttbar[14]*armttbar[45];
   armttbar[58]=armttbar[34]*armttbar[39];
   armttbar[45]=armttbar[52] + armttbar[58] + armttbar[45];
   armttbar[52]= - armttbar[63] + 1./2.*armttbar[67];
   armttbar[52]=armttbar[18]*armttbar[52];
   armttbar[58]=armttbar[27]*armttbar[55];
   armttbar[60]=armttbar[30]*armttbar[69];
   armttbar[42]=armttbar[25]*armttbar[42];
   armttbar[42]=armttbar[42] + armttbar[58] - 2./3.*armttbar[60];
   armttbar[42]=armttbar[42]*armttbar[68];
   armttbar[53]=MMt*armttbar[53];
   armttbar[53]=17*armttbar[63] + 5*armttbar[53];
   armttbar[58]= - armttbar[8]*armttbar[39];
   armttbar[53]=1./6.*armttbar[53] + armttbar[58];
   armttbar[53]=armttbar[22]*armttbar[53];
   armttbar[37]= - 14./27.*armttbar[37] - 112./27. + armttbar[39];
   armttbar[37]=armttbar[35]*armttbar[37];
   armttbar[39]=armttbar[33]*armttbar[57];
   armttbar[43]=armttbar[31]*armttbar[43]*armttbar[50];
   armttbar[37]=armttbar[44] + armttbar[59] + armttbar[53] + 
   armttbar[47] + armttbar[56] + armttbar[40] + armttbar[49] + 
   armttbar[43] + armttbar[48] + armttbar[38] + armttbar[41] + 
   armttbar[39] + armttbar[42] + 5./3.*armttbar[52] + armttbar[46] + 2*
   armttbar[37] + 2./3.*armttbar[45] + 1./3.*armttbar[51];
   armttbar[37]=armttbar[7]*armttbar[37];
   armttbar[38]=MMt*armttbar[4];
   armttbar[38]=armttbar[65] + 8*armttbar[38];
   armttbar[38]=MMt*armttbar[38];
   armttbar[38]= - 9*armttbar[54] + armttbar[38];
   armttbar[38]=armttbar[3]*armttbar[38]*armttbar[1]*armttbar[55];

      mttbarret = armttbar[37] + 4*armttbar[38];
      return mttbarret;
}
