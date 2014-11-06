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

#include <HH.hpp>
std::complex<long double>
HH<OS>::ygl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHHGL[58], yuHHGLret;

    aryuHHGL[1]=double(boson);
    aryuHHGL[2]=pow(SW,-1);
    aryuHHGL[3]=pow(MMH,-1);
    aryuHHGL[4]=pow(MMW,-1);
    aryuHHGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuHHGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuHHGL[7]=Tsil::I2(0,MMH,MMt,mu2);
    aryuHHGL[8]=Tsil::I2(0,0,MMH,mu2);
    aryuHHGL[9]=Tsil::I2(0,0,MMt,mu2);
    aryuHHGL[10]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHHGL[11]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHHGL[12]=Tsil::A(MMH,mu2);
    aryuHHGL[13]=Tsil::A(MMt,mu2);
    aryuHHGL[14]=std::real(Tsil::B(0,0,MMH,mu2));
    aryuHHGL[15]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuHHGL[16]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuHHGL[17]=Tsil::B(0,MMt,MMH,mu2);
    aryuHHGL[18]=Tsil::B(0,0,MMH,mu2);
    aryuHHGL[19]=Tsil::Beps(MMH,MMH,MMH,mu2);
    aryuHHGL[20]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHHGL[21]=pow(MMt,-1);
    aryuHHGL[22]=Tsil::Aeps(MMH,mu2);
    aryuHHGL[23]=Tsil::Aeps(MMt,mu2);
    aryuHHGL[24]=prottttt0->Suxv(0);
    aryuHHGL[25]=protHHHHH->M(0);
    aryuHHGL[26]=protHtHtt->M(0);
    aryuHHGL[27]=protttttH->M(0);
    aryuHHGL[28]=protHHHHH->Uzxyv(0);
    aryuHHGL[29]=protHtHtt->Uzxyv(0);
    aryuHHGL[30]=protHtHtt->Uuyxv(0);
    aryuHHGL[31]=protHtHtt->Svyz(0);
    aryuHHGL[32]=prot0H0H0->M(0);
    aryuHHGL[33]=prot0t0tt->M(0);
    aryuHHGL[34]=prot0t0t0->M(0);
    aryuHHGL[35]=prot0000H->M(0);
    aryuHHGL[36]=prot0H0H0->Uuyxv(0);
    aryuHHGL[37]=prot0t0t0->Uuyxv(0);
    aryuHHGL[38]=prot0H0H0->Tyzv(0);
    aryuHHGL[39]=prot0t0t0->Tyzv(0);
    aryuHHGL[40]=1/(4*MMt - MMH);
   aryuHHGL[41]=aryuHHGL[19] - aryuHHGL[29];
   aryuHHGL[42]=3./8. + aryuHHGL[11];
   aryuHHGL[42]=aryuHHGL[10]*aryuHHGL[42];
   aryuHHGL[42]=aryuHHGL[41] + aryuHHGL[42];
   aryuHHGL[43]=pow(Pi,2);
   aryuHHGL[44]=MMt*aryuHHGL[33];
   aryuHHGL[44]= - aryuHHGL[44] + 3./8.*aryuHHGL[43] + 3./4.*
   aryuHHGL[15] + 51./8. + aryuHHGL[39];
   aryuHHGL[45]=1./2.*aryuHHGL[11];
   aryuHHGL[46]= - 1./8.*aryuHHGL[11] - 1 - 1./2.*aryuHHGL[15];
   aryuHHGL[46]=aryuHHGL[46]*aryuHHGL[45];
   aryuHHGL[47]=aryuHHGL[17] + aryuHHGL[30];
   aryuHHGL[48]=3*aryuHHGL[11];
   aryuHHGL[49]= - 1 + aryuHHGL[48];
   aryuHHGL[49]=aryuHHGL[14]*aryuHHGL[49];
   aryuHHGL[42]=1./4.*aryuHHGL[49] + 3./8.*aryuHHGL[18] + aryuHHGL[46]
    + 1./2.*aryuHHGL[44] + 3*aryuHHGL[42] - 3./4.*aryuHHGL[47];
   aryuHHGL[42]=MMt*aryuHHGL[42];
   aryuHHGL[44]=3./4.*aryuHHGL[21] + aryuHHGL[40];
   aryuHHGL[44]=aryuHHGL[13]*aryuHHGL[44];
   aryuHHGL[44]=aryuHHGL[44] - 9*aryuHHGL[10] - aryuHHGL[45] - 7./2. + 
   aryuHHGL[15];
   aryuHHGL[44]=aryuHHGL[13]*aryuHHGL[44];
   aryuHHGL[46]=pow(MMt,2);
   aryuHHGL[47]=aryuHHGL[26]*aryuHHGL[46];
   aryuHHGL[49]=27./2.*aryuHHGL[10];
   aryuHHGL[50]=aryuHHGL[49] - 3 + 7./2.*aryuHHGL[14];
   aryuHHGL[50]=aryuHHGL[12]*aryuHHGL[50];
   aryuHHGL[42]=1./4.*aryuHHGL[44] + 1./16.*aryuHHGL[8] + 1./8.*
   aryuHHGL[50] + 9./2.*aryuHHGL[47] - 3./8.*aryuHHGL[7] + aryuHHGL[42]
   ;
   aryuHHGL[44]=aryuHHGL[45] - 1;
   aryuHHGL[45]=3./2.*aryuHHGL[11];
   aryuHHGL[47]=aryuHHGL[44]*aryuHHGL[45];
   aryuHHGL[50]= - 279./2. + 37./3.*aryuHHGL[43];
   aryuHHGL[47]=aryuHHGL[47] + 33*aryuHHGL[38] + 243*S2 + 1./2.*
   aryuHHGL[50] - 27*aryuHHGL[36];
   aryuHHGL[50]=3*aryuHHGL[13];
   aryuHHGL[51]= - 1 + aryuHHGL[11];
   aryuHHGL[51]=aryuHHGL[40]*aryuHHGL[51]*aryuHHGL[50];
   aryuHHGL[52]=7./4.*aryuHHGL[10] + 3./4. + aryuHHGL[14];
   aryuHHGL[52]=aryuHHGL[10]*aryuHHGL[52];
   aryuHHGL[52]=aryuHHGL[52] + aryuHHGL[19];
   aryuHHGL[53]= - 11 + aryuHHGL[14];
   aryuHHGL[53]=aryuHHGL[14]*aryuHHGL[53];
   aryuHHGL[44]=aryuHHGL[11]*aryuHHGL[44];
   aryuHHGL[44]=1./2. + aryuHHGL[44];
   aryuHHGL[44]=aryuHHGL[40]*aryuHHGL[44];
   aryuHHGL[44]=1./4.*aryuHHGL[44] + 27./2.*aryuHHGL[25] + 3*
   aryuHHGL[32] + 1./2.*aryuHHGL[35];
   aryuHHGL[44]=MMH*aryuHHGL[44];
   aryuHHGL[44]=3*aryuHHGL[44] + aryuHHGL[51] + 3./4.*aryuHHGL[53] + 9./
   2.*aryuHHGL[18] - 27./2.*aryuHHGL[28] + 1./2.*aryuHHGL[47] + 27*
   aryuHHGL[52];
   aryuHHGL[47]=1./8.*MMH;
   aryuHHGL[44]=aryuHHGL[44]*aryuHHGL[47];
   aryuHHGL[42]=3*aryuHHGL[42] + aryuHHGL[44];
   aryuHHGL[42]=9./32.*aryuHHGL[5] + 1./2.*aryuHHGL[42];
   aryuHHGL[42]=MMH*aryuHHGL[42];
   aryuHHGL[44]=1./4. - aryuHHGL[15];
   aryuHHGL[44]=5*aryuHHGL[44] - 7./8.*aryuHHGL[11];
   aryuHHGL[44]=aryuHHGL[44]*aryuHHGL[48];
   aryuHHGL[48]=19./2. - 3*aryuHHGL[39];
   aryuHHGL[48]=1./8.*aryuHHGL[48] + aryuHHGL[15];
   aryuHHGL[44]=aryuHHGL[44] + 6*aryuHHGL[30] + 3*aryuHHGL[48] - 29./32.
   *aryuHHGL[43];
   aryuHHGL[44]=aryuHHGL[3]*aryuHHGL[44];
   aryuHHGL[48]=3./2.*MMt;
   aryuHHGL[51]=pow(aryuHHGL[3],2);
   aryuHHGL[52]= - 9 - aryuHHGL[39];
   aryuHHGL[52]=aryuHHGL[52]*aryuHHGL[51]*aryuHHGL[48];
   aryuHHGL[44]=aryuHHGL[44] + aryuHHGL[52];
   aryuHHGL[44]=MMt*aryuHHGL[44];
   aryuHHGL[52]= - 1./2. - aryuHHGL[15];
   aryuHHGL[53]=7*aryuHHGL[15] - 1./4.*aryuHHGL[11];
   aryuHHGL[53]=aryuHHGL[11]*aryuHHGL[53];
   aryuHHGL[54]=aryuHHGL[3]*aryuHHGL[7];
   aryuHHGL[41]=aryuHHGL[44] - 9./8.*aryuHHGL[54] + 3./4.*aryuHHGL[53]
    + 3*aryuHHGL[52] - 35./64.*aryuHHGL[43] - 27./4.*aryuHHGL[41];
   aryuHHGL[41]=MMt*aryuHHGL[41];
   aryuHHGL[41]=27./16.*aryuHHGL[7] + aryuHHGL[41];
   aryuHHGL[41]=MMt*aryuHHGL[41];
   aryuHHGL[43]=pow(MMt,3);
   aryuHHGL[44]=2*MMt;
   aryuHHGL[52]= - aryuHHGL[3]*aryuHHGL[44];
   aryuHHGL[52]=1 + aryuHHGL[52];
   aryuHHGL[52]=aryuHHGL[52]*aryuHHGL[43];
   aryuHHGL[47]=aryuHHGL[46]*aryuHHGL[47];
   aryuHHGL[47]=aryuHHGL[52] + aryuHHGL[47];
   aryuHHGL[47]=aryuHHGL[27]*aryuHHGL[47];
   aryuHHGL[52]=aryuHHGL[44]*aryuHHGL[51];
   aryuHHGL[53]= - 2 - aryuHHGL[11];
   aryuHHGL[53]=aryuHHGL[53]*aryuHHGL[52];
   aryuHHGL[54]= - 5./8. + aryuHHGL[11];
   aryuHHGL[54]=aryuHHGL[3]*aryuHHGL[54];
   aryuHHGL[53]=aryuHHGL[54] + aryuHHGL[53];
   aryuHHGL[53]=MMt*aryuHHGL[53];
   aryuHHGL[54]= - 1 - aryuHHGL[11];
   aryuHHGL[53]=1./8.*aryuHHGL[54] + aryuHHGL[53];
   aryuHHGL[53]=MMt*aryuHHGL[53];
   aryuHHGL[51]=aryuHHGL[51]*MMt;
   aryuHHGL[54]=5./2.*aryuHHGL[3] - 7*aryuHHGL[51];
   aryuHHGL[54]=MMt*aryuHHGL[54];
   aryuHHGL[54]=11./4. + aryuHHGL[54];
   aryuHHGL[54]=aryuHHGL[12]*aryuHHGL[54];
   aryuHHGL[53]=aryuHHGL[53] + 3./16.*aryuHHGL[54];
   aryuHHGL[53]=aryuHHGL[12]*aryuHHGL[53];
   aryuHHGL[54]=5./8.*aryuHHGL[3];
   aryuHHGL[55]= - aryuHHGL[54] + aryuHHGL[52];
   aryuHHGL[55]=MMt*aryuHHGL[55];
   aryuHHGL[55]=3./4. + aryuHHGL[55];
   aryuHHGL[55]=MMt*aryuHHGL[55];
   aryuHHGL[55]=aryuHHGL[55] + 1./32.*MMH;
   aryuHHGL[55]=aryuHHGL[22]*aryuHHGL[55];
   aryuHHGL[56]=1./2.*aryuHHGL[3] - aryuHHGL[52];
   aryuHHGL[56]=MMt*aryuHHGL[56];
   aryuHHGL[56]= - 3./8. + aryuHHGL[56];
   aryuHHGL[56]=MMt*aryuHHGL[56];
   aryuHHGL[56]=aryuHHGL[56] + 3./16.*MMH;
   aryuHHGL[56]=aryuHHGL[6]*aryuHHGL[56];
   aryuHHGL[57]= - 9./4.*aryuHHGL[3] + 4*aryuHHGL[51];
   aryuHHGL[57]=MMt*aryuHHGL[57];
   aryuHHGL[57]=5./2. + aryuHHGL[57];
   aryuHHGL[57]=MMt*aryuHHGL[57];
   aryuHHGL[57]=aryuHHGL[57] + 1./4.*MMH;
   aryuHHGL[57]=aryuHHGL[23]*aryuHHGL[57];
   aryuHHGL[47]=aryuHHGL[56] + aryuHHGL[53] + aryuHHGL[55] + 
   aryuHHGL[57] + aryuHHGL[47];
   aryuHHGL[53]=aryuHHGL[54] - aryuHHGL[51];
   aryuHHGL[53]=MMt*aryuHHGL[53];
   aryuHHGL[53]=9./8. + aryuHHGL[53];
   aryuHHGL[53]=aryuHHGL[17]*aryuHHGL[53];
   aryuHHGL[54]=7./4.*aryuHHGL[3] - aryuHHGL[51];
   aryuHHGL[54]=aryuHHGL[9]*aryuHHGL[54];
   aryuHHGL[53]=aryuHHGL[54] + aryuHHGL[53];
   aryuHHGL[54]= - 1./2. - aryuHHGL[11];
   aryuHHGL[49]=aryuHHGL[54]*aryuHHGL[49];
   aryuHHGL[54]= - aryuHHGL[14]*aryuHHGL[45];
   aryuHHGL[55]=MMt*aryuHHGL[3];
   aryuHHGL[56]= - 1./2. + aryuHHGL[55];
   aryuHHGL[56]=aryuHHGL[37]*aryuHHGL[56];
   aryuHHGL[57]=aryuHHGL[31]*aryuHHGL[3];
   aryuHHGL[49]= - 51./8.*aryuHHGL[57] + 15./4.*aryuHHGL[56] + 
   aryuHHGL[49] + aryuHHGL[54] + 3./2.*aryuHHGL[53];
   aryuHHGL[49]=aryuHHGL[46]*aryuHHGL[49];
   aryuHHGL[53]=13./4.*aryuHHGL[11] + 119./16. - aryuHHGL[15];
   aryuHHGL[53]=aryuHHGL[3]*aryuHHGL[53];
   aryuHHGL[51]=aryuHHGL[53] + 9./2.*aryuHHGL[51];
   aryuHHGL[51]=MMt*aryuHHGL[51];
   aryuHHGL[45]= - aryuHHGL[45] + 57./4. - aryuHHGL[15];
   aryuHHGL[45]= - 9./4.*aryuHHGL[10] - 1./2.*aryuHHGL[14] + 1./4.*
   aryuHHGL[45] + aryuHHGL[51];
   aryuHHGL[45]=MMt*aryuHHGL[45];
   aryuHHGL[51]= - aryuHHGL[3] - aryuHHGL[52];
   aryuHHGL[44]=aryuHHGL[51]*aryuHHGL[44];
   aryuHHGL[44]=5./16. + aryuHHGL[44];
   aryuHHGL[44]=aryuHHGL[12]*aryuHHGL[44];
   aryuHHGL[51]= - 11./4.*aryuHHGL[3] - aryuHHGL[52];
   aryuHHGL[51]=MMt*aryuHHGL[51];
   aryuHHGL[51]=9./16. + aryuHHGL[51];
   aryuHHGL[51]=aryuHHGL[13]*aryuHHGL[51];
   aryuHHGL[44]=aryuHHGL[51] + aryuHHGL[44] + aryuHHGL[45];
   aryuHHGL[44]=aryuHHGL[44]*aryuHHGL[50];
   aryuHHGL[45]=5./2. - 13*aryuHHGL[55];
   aryuHHGL[45]=aryuHHGL[45]*aryuHHGL[46];
   aryuHHGL[50]=MMH*aryuHHGL[48];
   aryuHHGL[45]=aryuHHGL[45] + aryuHHGL[50];
   aryuHHGL[45]=aryuHHGL[20]*aryuHHGL[45];
   aryuHHGL[50]=1 - 5*aryuHHGL[11];
   aryuHHGL[50]=aryuHHGL[50]*aryuHHGL[55];
   aryuHHGL[51]= - 3./2. + aryuHHGL[11];
   aryuHHGL[50]=1./2.*aryuHHGL[51] + aryuHHGL[50];
   aryuHHGL[46]=aryuHHGL[50]*aryuHHGL[46];
   aryuHHGL[50]= - 1./2. - aryuHHGL[55];
   aryuHHGL[50]=aryuHHGL[13]*MMt*aryuHHGL[50];
   aryuHHGL[46]=aryuHHGL[46] + aryuHHGL[50];
   aryuHHGL[46]=aryuHHGL[16]*aryuHHGL[46];
   aryuHHGL[45]=aryuHHGL[45] + aryuHHGL[46];
   aryuHHGL[46]= - 3./2.*aryuHHGL[34] - 18*aryuHHGL[26];
   aryuHHGL[43]=aryuHHGL[43]*aryuHHGL[46];
   aryuHHGL[46]= - 1 - 1./4.*aryuHHGL[55];
   aryuHHGL[46]=aryuHHGL[24]*aryuHHGL[46]*aryuHHGL[48];
   aryuHHGL[41]=aryuHHGL[46] + aryuHHGL[44] + aryuHHGL[41] + 
   aryuHHGL[43] + aryuHHGL[49] + aryuHHGL[42] + 3./4.*aryuHHGL[45] + 3*
   aryuHHGL[47];

      yuHHGLret = aryuHHGL[41]*pow(aryuHHGL[4],2)*pow(aryuHHGL[2],4)*
      aryuHHGL[1];
      return yuHHGLret;
}
