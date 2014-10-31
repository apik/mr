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
tt::mgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[51], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=pow(SW,-1);
    armttbarGL[3]=pow(MMH,-1);
    armttbarGL[4]=pow(MMW,-1);
    armttbarGL[5]=pow(MMt,-1);
    armttbarGL[6]=Tsil::I2(MMH,MMH,MMH,mu2);
    armttbarGL[7]=Tsil::I2(MMH,MMt,MMt,mu2);
    armttbarGL[8]=Tsil::I2(0,0,MMH,mu2);
    armttbarGL[9]=Tsil::I2(0,0,MMt,mu2);
    armttbarGL[10]=Tsil::B(MMH,MMH,MMH,mu2);
    armttbarGL[11]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbarGL[12]=Tsil::A(MMH,mu2);
    armttbarGL[13]=Tsil::A(MMt,mu2);
    armttbarGL[14]=Tsil::B(MMt,MMt,MMH,mu2);
    armttbarGL[15]=std::real(Tsil::B(0,0,MMH,mu2));
    armttbarGL[16]=std::real(Tsil::B(0,0,MMt,mu2));
    armttbarGL[17]=Tsil::B(0,MMH,MMt,mu2);
    armttbarGL[18]=Tsil::Beps(MMH,MMt,MMt,mu2);
    armttbarGL[19]=Tsil::Aeps(MMH,mu2);
    armttbarGL[20]=Tsil::Aeps(MMt,mu2);
    armttbarGL[21]=prot0ttHt->Tuxv(0);
    armttbarGL[22]=protHHttH->M(0);
    armttbarGL[23]=protHttHt->M(0);
    armttbarGL[24]=protHHttH->Uxzuv(0);
    armttbarGL[25]=protHHttH->Suxv(0);
    armttbarGL[26]=protHHttH->Uuyxv(0);
    armttbarGL[27]=protWt000->Tyzv(0);
    armttbarGL[28]=protH0tt0->M(0);
    armttbarGL[29]=protH0t00->M(0);
    armttbarGL[30]=prot0Htt0->M(0);
    armttbarGL[31]=prot0H0t0->M(0);
    armttbarGL[32]=prot00ttH->M(0);
    armttbarGL[33]=protH0t00->Uxzuv(0);
    armttbarGL[34]=protH0t00->Uzxyv(0);
    armttbarGL[35]=protH0tt0->Uuyxv(0);
    armttbarGL[36]=protH0t00->Txuv(0);
   armttbarGL[37]=3*armttbarGL[11];
   armttbarGL[38]= - 19./4. + armttbarGL[11];
   armttbarGL[38]=armttbarGL[38]*armttbarGL[37];
   armttbarGL[39]=pow(Pi,2);
   armttbarGL[39]=armttbarGL[39] - armttbarGL[17];
   armttbarGL[40]=armttbarGL[29] + armttbarGL[31];
   armttbarGL[40]=3*armttbarGL[22] - armttbarGL[23] + 1./4.*
   armttbarGL[40];
   armttbarGL[40]=MMH*armttbarGL[40];
   armttbarGL[38]=armttbarGL[40] + 75./16.*armttbarGL[27] - 115./8.*
   armttbarGL[21] + armttbarGL[38] + 1./4.*armttbarGL[36] - 473./128.
    - armttbarGL[24] + 3./32.*armttbarGL[39];
   armttbarGL[39]=3./2.*armttbarGL[3];
   armttbarGL[40]=pow(armttbarGL[13],2);
   armttbarGL[39]= - armttbarGL[40]*armttbarGL[39];
   armttbarGL[41]=armttbarGL[11]*armttbarGL[13];
   armttbarGL[42]=3./2.*armttbarGL[7] + armttbarGL[9];
   armttbarGL[39]=armttbarGL[39] + 1./2.*armttbarGL[42] - 3*
   armttbarGL[41];
   armttbarGL[42]=3*armttbarGL[3];
   armttbarGL[39]=armttbarGL[39]*armttbarGL[42];
   armttbarGL[43]=armttbarGL[12]*armttbarGL[3];
   armttbarGL[44]=armttbarGL[3]*armttbarGL[13];
   armttbarGL[45]= - 7 + 5*armttbarGL[11];
   armttbarGL[45]= - 1./4.*armttbarGL[43] + 1./8.*armttbarGL[45] + 
   armttbarGL[44];
   armttbarGL[45]=armttbarGL[14]*armttbarGL[45];
   armttbarGL[37]=5./2. - armttbarGL[37];
   armttbarGL[37]=armttbarGL[37]*armttbarGL[42];
   armttbarGL[46]=3./2.*armttbarGL[16];
   armttbarGL[47]= - armttbarGL[3]*armttbarGL[46];
   armttbarGL[37]=armttbarGL[47] + armttbarGL[23] + armttbarGL[37];
   armttbarGL[47]=armttbarGL[11] - 1;
   armttbarGL[48]= - 1./2.*armttbarGL[47] - 6*armttbarGL[44];
   armttbarGL[42]=armttbarGL[14]*armttbarGL[48]*armttbarGL[42];
   armttbarGL[37]=1./2.*armttbarGL[37] + armttbarGL[42];
   armttbarGL[37]=MMt*armttbarGL[37];
   armttbarGL[42]=1./2.*armttbarGL[12];
   armttbarGL[48]=1./8.*armttbarGL[16];
   armttbarGL[49]=armttbarGL[48] - 33./8. + armttbarGL[11];
   armttbarGL[49]=armttbarGL[42]*armttbarGL[49];
   armttbarGL[49]= - 6*armttbarGL[20] - 185./32.*armttbarGL[19] + 
   armttbarGL[49];
   armttbarGL[49]=armttbarGL[3]*armttbarGL[49];
   armttbarGL[50]=3./4. + 11*armttbarGL[11];
   armttbarGL[50]=3./16.*armttbarGL[16] + 1./4.*armttbarGL[50] - 15*
   armttbarGL[44];
   armttbarGL[48]=armttbarGL[50]*armttbarGL[48];
   armttbarGL[37]=armttbarGL[37] - 5./32.*armttbarGL[33] + 3*
   armttbarGL[45] + 13./32.*armttbarGL[18] + armttbarGL[48] + 1./4.*
   armttbarGL[38] + armttbarGL[39] + armttbarGL[49];
   armttbarGL[38]=pow(armttbarGL[2],4)*pow(armttbarGL[4],2);
   armttbarGL[37]=MMt*armttbarGL[38]*armttbarGL[37];
   armttbarGL[39]=9*armttbarGL[10];
   armttbarGL[45]=3*armttbarGL[36] + 89 - 3./2.*armttbarGL[17];
   armttbarGL[48]=armttbarGL[28] + armttbarGL[30] + armttbarGL[32] - 9*
   armttbarGL[22] - armttbarGL[23];
   armttbarGL[48]=MMH*armttbarGL[48];
   armttbarGL[45]=armttbarGL[48] - armttbarGL[39] + 1./2.*
   armttbarGL[45] - 9*armttbarGL[26];
   armttbarGL[48]=211./2. + armttbarGL[39];
   armttbarGL[48]=1./8.*armttbarGL[48] - armttbarGL[11];
   armttbarGL[48]=armttbarGL[11]*armttbarGL[48];
   armttbarGL[47]=armttbarGL[15]*armttbarGL[47];
   armttbarGL[45]= - 11./8.*armttbarGL[34] + 5./8.*armttbarGL[27] + 3./
   8.*armttbarGL[47] + 35./4.*armttbarGL[21] + armttbarGL[48] + 1./8.*
   armttbarGL[45] + 5./16.*armttbarGL[33] + 35./16.*armttbarGL[18];
   armttbarGL[45]=MMH*armttbarGL[45];
   armttbarGL[47]=9*armttbarGL[15] + 27*armttbarGL[10];
   armttbarGL[47]=armttbarGL[13]*armttbarGL[47];
   armttbarGL[48]= - 145./16.*armttbarGL[13] + armttbarGL[7];
   armttbarGL[47]=33*armttbarGL[41] + 3*armttbarGL[48] + 5./4.*
   armttbarGL[9] + armttbarGL[47];
   armttbarGL[44]=armttbarGL[46] - 207./4. + 203*armttbarGL[44];
   armttbarGL[43]=1./8.*armttbarGL[44] + armttbarGL[43];
   armttbarGL[43]=armttbarGL[12]*armttbarGL[43];
   armttbarGL[44]=armttbarGL[11] - 3./2.;
   armttbarGL[46]= - MMH*armttbarGL[44];
   armttbarGL[42]=armttbarGL[42] + armttbarGL[13] + armttbarGL[46];
   armttbarGL[42]=armttbarGL[14]*armttbarGL[42];
   armttbarGL[46]=armttbarGL[16]*armttbarGL[13];
   armttbarGL[42]=3./2.*armttbarGL[42] + 17./16.*armttbarGL[25] + 85./
   16.*armttbarGL[19] + armttbarGL[43] + 3./4.*armttbarGL[46] + 1./2.*
   armttbarGL[47] + armttbarGL[45];
   armttbarGL[42]=armttbarGL[42]*armttbarGL[38];
   armttbarGL[37]=1./4.*armttbarGL[42] + armttbarGL[37];
   armttbarGL[37]=MMt*armttbarGL[37];
   armttbarGL[42]=MMH*armttbarGL[5];
   armttbarGL[43]= - armttbarGL[36] - 1 - armttbarGL[17];
   armttbarGL[43]=armttbarGL[43]*armttbarGL[42];
   armttbarGL[45]=armttbarGL[5]*armttbarGL[13];
   armttbarGL[43]=armttbarGL[35] + armttbarGL[43] + 5./2.*
   armttbarGL[17] + 3./2.*armttbarGL[24] - 39./2. + armttbarGL[45];
   armttbarGL[45]=1./2.*armttbarGL[11];
   armttbarGL[46]= - 5 - 3*armttbarGL[10];
   armttbarGL[46]=3*armttbarGL[46] + armttbarGL[11];
   armttbarGL[46]=armttbarGL[46]*armttbarGL[45];
   armttbarGL[47]=1./2.*armttbarGL[5];
   armttbarGL[48]= - armttbarGL[8]*armttbarGL[47];
   armttbarGL[44]=armttbarGL[15]*armttbarGL[44];
   armttbarGL[43]= - 3./2.*armttbarGL[44] + armttbarGL[48] - 23./4.*
   armttbarGL[21] + armttbarGL[46] + 27./4.*armttbarGL[10] + 3*
   armttbarGL[26] + 1./2.*armttbarGL[43];
   armttbarGL[43]=MMH*armttbarGL[43];
   armttbarGL[39]=3*armttbarGL[15] + armttbarGL[39] + 31;
   armttbarGL[39]=armttbarGL[13]*armttbarGL[39];
   armttbarGL[44]=armttbarGL[5]*armttbarGL[40];
   armttbarGL[39]= - 5./2.*armttbarGL[8] - 23*armttbarGL[41] - 13*
   armttbarGL[7] - 3*armttbarGL[44] + armttbarGL[39];
   armttbarGL[39]=1./2.*armttbarGL[39] + armttbarGL[43];
   armttbarGL[39]=MMH*armttbarGL[39];
   armttbarGL[41]=MMH - armttbarGL[13];
   armttbarGL[41]=armttbarGL[47]*armttbarGL[41];
   armttbarGL[41]=3./4.*armttbarGL[15] - armttbarGL[45] + 9./4.*
   armttbarGL[10] + 9 + armttbarGL[41];
   armttbarGL[41]=MMH*armttbarGL[41];
   armttbarGL[41]=1./8.*armttbarGL[12] - 51./2.*armttbarGL[13] + 
   armttbarGL[41];
   armttbarGL[41]=armttbarGL[12]*armttbarGL[41];
   armttbarGL[39]=armttbarGL[39] + armttbarGL[41];
   armttbarGL[41]=31 + armttbarGL[42];
   armttbarGL[41]=armttbarGL[20]*armttbarGL[41];
   armttbarGL[41]=17./64.*armttbarGL[19] + 1./32.*armttbarGL[41] - 3./
   64.*armttbarGL[6];
   armttbarGL[41]=MMH*armttbarGL[41];
   armttbarGL[42]= - 21./64.*armttbarGL[18] + 3./32.*armttbarGL[34];
   armttbarGL[42]=armttbarGL[42]*pow(MMH,2);
   armttbarGL[39]=armttbarGL[40] + 1./16.*armttbarGL[39] + 
   armttbarGL[42] + armttbarGL[41];
   armttbarGL[38]=armttbarGL[39]*armttbarGL[38];
   armttbarGL[37]=armttbarGL[37] + armttbarGL[38];

      mttbarGLret = armttbarGL[37]*armttbarGL[1];
      return mttbarGLret;
}
