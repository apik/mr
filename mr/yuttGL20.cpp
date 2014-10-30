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
tt::mygl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[54], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=pow(SW,-1);
    aryuttGL[3]=pow(MMH,-1);
    aryuttGL[4]=pow(MMW,-1);
    aryuttGL[5]=pow(MMt,-1);
    aryuttGL[6]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuttGL[7]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuttGL[8]=Tsil::I2(0,MMH,MMt,mu2);
    aryuttGL[9]=Tsil::I2(0,0,MMH,mu2);
    aryuttGL[10]=Tsil::I2(0,0,MMt,mu2);
    aryuttGL[11]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuttGL[12]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[13]=Tsil::A(MMH,mu2);
    aryuttGL[14]=Tsil::A(MMt,mu2);
    aryuttGL[15]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuttGL[16]=std::real(Tsil::B(0,0,MMH,mu2));
    aryuttGL[17]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuttGL[18]=Tsil::B(0,MMH,MMt,mu2);
    aryuttGL[19]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryuttGL[20]=Tsil::Aeps(MMH,mu2);
    aryuttGL[21]=Tsil::Aeps(MMt,mu2);
    aryuttGL[22]=prot0ttHt->Tuxv(0);
    aryuttGL[23]=protHHttH->M(0);
    aryuttGL[24]=protHttHt->M(0);
    aryuttGL[25]=protHHttH->Uxzuv(0);
    aryuttGL[26]=protHHttH->Suxv(0);
    aryuttGL[27]=protHHttH->Uuyxv(0);
    aryuttGL[28]=protWt000->Tyzv(0);
    aryuttGL[29]=protH0tt0->M(0);
    aryuttGL[30]=protH0t00->M(0);
    aryuttGL[31]=prot0Htt0->M(0);
    aryuttGL[32]=prot0H0t0->M(0);
    aryuttGL[33]=prot00ttH->M(0);
    aryuttGL[34]=protH0t00->Uxzuv(0);
    aryuttGL[35]=protH0t00->Uzxyv(0);
    aryuttGL[36]=protH0tt0->Uuyxv(0);
    aryuttGL[37]=protH0t00->Txuv(0);
   aryuttGL[38]=3*aryuttGL[12];
   aryuttGL[39]= - 15./2. + aryuttGL[12];
   aryuttGL[39]=aryuttGL[39]*aryuttGL[38];
   aryuttGL[40]=1./8.*aryuttGL[17];
   aryuttGL[41]=3./4.*aryuttGL[17] - 39./4. + 11*aryuttGL[12];
   aryuttGL[41]=aryuttGL[41]*aryuttGL[40];
   aryuttGL[39]=1./4.*aryuttGL[37] + 75./16.*aryuttGL[28] + 
   aryuttGL[41] - 1073./128. + aryuttGL[39];
   aryuttGL[41]=pow(aryuttGL[13],2);
   aryuttGL[42]=aryuttGL[41]*aryuttGL[3];
   aryuttGL[43]=aryuttGL[10] - 3./8.*aryuttGL[8];
   aryuttGL[40]=aryuttGL[40] - 3./8. + aryuttGL[12];
   aryuttGL[40]=aryuttGL[13]*aryuttGL[40];
   aryuttGL[40]= - 63./16.*aryuttGL[42] + aryuttGL[40] - 9./4.*
   aryuttGL[7] + 3*aryuttGL[43] - 185./16.*aryuttGL[20];
   aryuttGL[40]=aryuttGL[3]*aryuttGL[40];
   aryuttGL[43]=aryuttGL[3]*aryuttGL[13];
   aryuttGL[44]=aryuttGL[14]*aryuttGL[3];
   aryuttGL[45]= - 1 + 5./2.*aryuttGL[12];
   aryuttGL[45]= - aryuttGL[44] + 1./2.*aryuttGL[45] + aryuttGL[43];
   aryuttGL[46]=3*aryuttGL[15];
   aryuttGL[45]=aryuttGL[45]*aryuttGL[46];
   aryuttGL[47]=7./8. - aryuttGL[12];
   aryuttGL[47]=aryuttGL[47]*aryuttGL[44];
   aryuttGL[39]= - 115./16.*aryuttGL[22] + aryuttGL[45] + 3*
   aryuttGL[47] + 1./2.*aryuttGL[39] + aryuttGL[40];
   aryuttGL[40]=aryuttGL[30] + aryuttGL[32];
   aryuttGL[45]=1./16.*MMH;
   aryuttGL[40]=aryuttGL[45]*aryuttGL[40];
   aryuttGL[45]=aryuttGL[12] - 1;
   aryuttGL[47]= - aryuttGL[15]*aryuttGL[45];
   aryuttGL[47]=aryuttGL[47] - 37./16. - aryuttGL[12];
   aryuttGL[48]=pow(Pi,2);
   aryuttGL[47]= - 29./32.*aryuttGL[48] + 3*aryuttGL[47];
   aryuttGL[47]=aryuttGL[3]*aryuttGL[47];
   aryuttGL[47]=aryuttGL[24] + aryuttGL[47];
   aryuttGL[47]=MMt*aryuttGL[47];
   aryuttGL[49]=aryuttGL[48] + aryuttGL[25];
   aryuttGL[50]=aryuttGL[21]*aryuttGL[3];
   aryuttGL[51]=1./4.*MMH;
   aryuttGL[52]= - aryuttGL[24]*aryuttGL[51];
   aryuttGL[53]=aryuttGL[23]*MMH;
   aryuttGL[39]=1./2.*aryuttGL[47] + 13./32.*aryuttGL[19] - 5./32.*
   aryuttGL[34] - 3./128.*aryuttGL[18] + 3./4.*aryuttGL[53] + 
   aryuttGL[52] + 1./2.*aryuttGL[39] - 6*aryuttGL[50] - 1./4.*
   aryuttGL[49] + aryuttGL[40];
   aryuttGL[40]=pow(aryuttGL[2],4)*pow(aryuttGL[4],2);
   aryuttGL[39]=MMt*aryuttGL[40]*aryuttGL[39];
   aryuttGL[47]=9./8.*aryuttGL[11] + 3./8.*aryuttGL[16];
   aryuttGL[45]=aryuttGL[45]*aryuttGL[47];
   aryuttGL[47]=117./8. - aryuttGL[12];
   aryuttGL[47]=aryuttGL[12]*aryuttGL[47];
   aryuttGL[49]=aryuttGL[12] - 1./4.;
   aryuttGL[50]=aryuttGL[15]*aryuttGL[49];
   aryuttGL[52]=MMH*aryuttGL[33];
   aryuttGL[45]= - 11./8.*aryuttGL[35] + 1./8.*aryuttGL[52] + 35./4.*
   aryuttGL[22] - 3./2.*aryuttGL[50] + 3./16.*aryuttGL[37] + 5./8.*
   aryuttGL[28] + 1./32.*aryuttGL[17] + 361./32. + aryuttGL[47] + 
   aryuttGL[45] + 9./16.*aryuttGL[48] + 35./16.*aryuttGL[19] - 9./8.*
   aryuttGL[27] + 5./16.*aryuttGL[34] - 3./32.*aryuttGL[18];
   aryuttGL[45]=MMH*aryuttGL[45];
   aryuttGL[47]=85./2.*aryuttGL[20] + 5*aryuttGL[10] + 27*aryuttGL[8];
   aryuttGL[50]= - 87./8. - aryuttGL[12];
   aryuttGL[50]=aryuttGL[13]*aryuttGL[50];
   aryuttGL[42]=61./4.*aryuttGL[42] + 3*aryuttGL[50] + 1./2.*
   aryuttGL[47] + 15*aryuttGL[7];
   aryuttGL[38]= - 1./4.*aryuttGL[17] - 133./16. + aryuttGL[38];
   aryuttGL[38]= - 93./2.*aryuttGL[44] + 3*aryuttGL[38] + 125./4.*
   aryuttGL[43];
   aryuttGL[38]=aryuttGL[14]*aryuttGL[38];
   aryuttGL[43]=aryuttGL[14] - aryuttGL[13];
   aryuttGL[44]=aryuttGL[43]*aryuttGL[46];
   aryuttGL[38]=aryuttGL[44] + 1./2.*aryuttGL[42] + aryuttGL[38];
   aryuttGL[42]=aryuttGL[29] - aryuttGL[24];
   aryuttGL[44]=pow(MMH,2);
   aryuttGL[46]=1./8.*aryuttGL[44];
   aryuttGL[42]=aryuttGL[46]*aryuttGL[42];
   aryuttGL[47]=aryuttGL[23]*aryuttGL[44];
   aryuttGL[46]=aryuttGL[31]*aryuttGL[46];
   aryuttGL[38]=aryuttGL[46] - 9./8.*aryuttGL[47] + 17./16.*
   aryuttGL[26] + 1./2.*aryuttGL[38] + aryuttGL[42] + aryuttGL[45];
   aryuttGL[38]=aryuttGL[40]*aryuttGL[38];
   aryuttGL[38]=1./4.*aryuttGL[38] + aryuttGL[39];
   aryuttGL[38]=MMt*aryuttGL[1]*aryuttGL[38];
   aryuttGL[39]=3*aryuttGL[16];
   aryuttGL[42]= - aryuttGL[49]*aryuttGL[39];
   aryuttGL[45]=MMH*aryuttGL[5];
   aryuttGL[46]= - 1 - aryuttGL[37];
   aryuttGL[46]=aryuttGL[46]*aryuttGL[45];
   aryuttGL[47]=aryuttGL[21] + aryuttGL[14] + aryuttGL[13];
   aryuttGL[47]=aryuttGL[5]*aryuttGL[47];
   aryuttGL[50]= - 61./4. + aryuttGL[12];
   aryuttGL[50]=aryuttGL[12]*aryuttGL[50];
   aryuttGL[42]=aryuttGL[46] - 23./2.*aryuttGL[22] + aryuttGL[36] + 
   aryuttGL[42] - 571./16. + aryuttGL[50] + aryuttGL[47];
   aryuttGL[42]=aryuttGL[42]*aryuttGL[51];
   aryuttGL[46]=131./2. + aryuttGL[12];
   aryuttGL[39]=1./2.*aryuttGL[46] - aryuttGL[39];
   aryuttGL[39]=aryuttGL[13]*aryuttGL[39];
   aryuttGL[39]=aryuttGL[39] - 9*aryuttGL[8] + 17./2.*aryuttGL[20];
   aryuttGL[46]= - MMH*aryuttGL[49];
   aryuttGL[43]=aryuttGL[46] + aryuttGL[43];
   aryuttGL[43]=aryuttGL[11]*aryuttGL[43];
   aryuttGL[46]=3./4.*aryuttGL[14] - 1./2.*aryuttGL[13];
   aryuttGL[46]=aryuttGL[5]*aryuttGL[46];
   aryuttGL[46]=3./2.*aryuttGL[16] + 27./2. - 7*aryuttGL[12] + 
   aryuttGL[46];
   aryuttGL[46]=aryuttGL[14]*aryuttGL[46];
   aryuttGL[39]=9./4.*aryuttGL[43] + aryuttGL[42] + 31./4.*aryuttGL[21]
    + 1./2.*aryuttGL[46] - 3./8.*aryuttGL[6] - aryuttGL[7] + 1./4.*
   aryuttGL[39];
   aryuttGL[39]=MMH*aryuttGL[39];
   aryuttGL[42]=5./2. - aryuttGL[45];
   aryuttGL[42]=aryuttGL[18]*aryuttGL[42];
   aryuttGL[42]=3./8.*aryuttGL[25] + 3./4.*aryuttGL[35] + 243./8.*S2 - 
   1./6.*aryuttGL[48] - 21./8.*aryuttGL[19] + 3./2.*aryuttGL[27] + 1./4.
   *aryuttGL[42];
   aryuttGL[42]=aryuttGL[44]*aryuttGL[42];
   aryuttGL[43]= - 9*aryuttGL[13] + 25*aryuttGL[14];
   aryuttGL[43]=aryuttGL[14]*aryuttGL[43];
   aryuttGL[41]= - 19./4.*aryuttGL[41] + aryuttGL[43];
   aryuttGL[43]= - 5./2. - aryuttGL[45];
   aryuttGL[43]=aryuttGL[9]*aryuttGL[43]*aryuttGL[51];
   aryuttGL[39]=aryuttGL[43] + 1./2.*aryuttGL[41] + aryuttGL[39] + 
   aryuttGL[42];
   aryuttGL[41]=1./8.*aryuttGL[1];
   aryuttGL[39]=aryuttGL[41]*aryuttGL[39]*aryuttGL[40];

      yuttGLret = aryuttGL[38] + aryuttGL[39];
      return yuttGLret;
}
