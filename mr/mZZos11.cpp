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
ZZ<MS>::x11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZos[34], mZZosret;

    armZZos[1]=double(nH);
    armZZos[2]=pow(mmZ,-1);
    armZZos[3]=pow(s,-1);
    armZZos[4]=pow(c,-1);
    armZZos[5]=pow(mmH,-1);
    armZZos[6]=Tsil::B(mmt,mmt,mmZ,mu2);
    armZZos[7]=Tsil::A(mmt,mu2);
    armZZos[8]=Tsil::Beps(mmt,mmt,mmZ,mu2);
    armZZos[9]=Tsil::Aeps(mmt,mu2);
    armZZos[10]=std::real(Tsil::B(0,0,mmZ,mu2));
    armZZos[11]=prottttt0->M(0);
    armZZos[12]=prot00000->M(0);
    armZZos[13]=prottttt0->Vzxyv(0);
    armZZos[14]=prottttt0->Suxv(0);
    armZZos[15]=double(nL);
    armZZos[16]=1/(4*mmt - mmZ);
   armZZos[17]=pow(armZZos[3],2);
   armZZos[18]=pow(armZZos[4],2);
   armZZos[19]=armZZos[17] + armZZos[18];
   armZZos[20]=armZZos[19]*pow(armZZos[2],2);
   armZZos[21]=4*mmt;
   armZZos[22]=armZZos[20]*armZZos[21];
   armZZos[23]=7./9.*armZZos[2];
   armZZos[23]=armZZos[23]*armZZos[18];
   armZZos[24]=armZZos[17] + 64./9.;
   armZZos[24]=armZZos[24]*armZZos[2];
   armZZos[23]=armZZos[24] - armZZos[23];
   armZZos[24]= - 5*armZZos[23] - armZZos[22];
   armZZos[25]=25./9.*armZZos[18] + armZZos[17] - 64./9.;
   armZZos[26]=armZZos[16]*armZZos[25];
   armZZos[24]=10*armZZos[26] + 2./3.*armZZos[24];
   armZZos[24]=armZZos[7]*armZZos[24];
   armZZos[26]= - 13./3.*armZZos[18] + 64./3. + armZZos[17];
   armZZos[26]=armZZos[2]*armZZos[26];
   armZZos[22]=5*armZZos[26] - armZZos[22];
   armZZos[26]=2./3.*mmt;
   armZZos[22]=armZZos[22]*armZZos[26];
   armZZos[27]=armZZos[20]*mmt;
   armZZos[27]= - 2*armZZos[27] + armZZos[23];
   armZZos[27]=armZZos[27]*armZZos[26];
   armZZos[28]=mmZ*armZZos[16];
   armZZos[29]=armZZos[28] + 1;
   armZZos[29]=armZZos[29]*armZZos[25];
   armZZos[27]=armZZos[27] + armZZos[29];
   armZZos[27]=armZZos[6]*armZZos[27];
   armZZos[30]=1./2.*armZZos[17] - 32./9. + 25./18.*armZZos[18];
   armZZos[28]=armZZos[30]*armZZos[28];
   armZZos[22]=armZZos[27] - 7*armZZos[28] + armZZos[22] - 311./18.*
   armZZos[18] + 352./9. - 15./2.*armZZos[17] + armZZos[24];
   armZZos[22]=armZZos[6]*armZZos[22];
   armZZos[24]=1./3.*armZZos[20];
   armZZos[27]=armZZos[19]*armZZos[2];
   armZZos[28]= - armZZos[5]*armZZos[27];
   armZZos[28]= - armZZos[24] + armZZos[28];
   armZZos[28]=armZZos[28]*armZZos[21];
   armZZos[31]= - 55./9.*armZZos[18] + 256./9. + armZZos[17];
   armZZos[31]=armZZos[2]*armZZos[31];
   armZZos[32]=12*armZZos[5];
   armZZos[27]=armZZos[27]*armZZos[32];
   armZZos[32]= - armZZos[7]*armZZos[27];
   armZZos[28]=armZZos[32] + 1./3.*armZZos[31] + armZZos[28];
   armZZos[28]=armZZos[7]*armZZos[28];
   armZZos[25]=armZZos[9]*armZZos[25];
   armZZos[31]=25./3.*armZZos[18] - 64./3. + 3*armZZos[17];
   armZZos[32]=armZZos[7]*armZZos[2];
   armZZos[32]=2*armZZos[32] - 1;
   armZZos[31]=armZZos[7]*armZZos[31]*armZZos[32];
   armZZos[25]=armZZos[25] + armZZos[31];
   armZZos[25]=armZZos[16]*armZZos[25];
   armZZos[31]=17./9.*armZZos[18] + armZZos[17] - 32./9.;
   armZZos[32]=armZZos[2]*armZZos[31];
   armZZos[33]=mmt*armZZos[24];
   armZZos[32]=armZZos[33] + armZZos[32];
   armZZos[32]=armZZos[14]*armZZos[32];
   armZZos[25]=armZZos[32] + armZZos[28] + armZZos[25];
   armZZos[28]=85./9.*armZZos[18] + 128./9. + 13*armZZos[17];
   armZZos[28]=armZZos[2]*armZZos[28];
   armZZos[20]=2*armZZos[20];
   armZZos[20]= - armZZos[9]*armZZos[20];
   armZZos[20]=armZZos[20] + armZZos[28];
   armZZos[24]=armZZos[24] - armZZos[27];
   armZZos[24]=mmt*armZZos[24];
   armZZos[26]= - armZZos[23]*armZZos[26];
   armZZos[19]=armZZos[26] + armZZos[19];
   armZZos[19]=armZZos[11]*armZZos[19];
   armZZos[19]=armZZos[19] + 1./3.*armZZos[20] + armZZos[24];
   armZZos[19]=armZZos[21]*armZZos[19];
   armZZos[20]= - 95./9.*armZZos[18] + 128./9. - 7*armZZos[17];
   armZZos[20]=armZZos[9]*armZZos[2]*armZZos[20];
   armZZos[20]=4*armZZos[20] - 431./9.*armZZos[18] + 764./9. - 37*
   armZZos[17];
   armZZos[21]=armZZos[17] - 8./9. + 5./9.*armZZos[18];
   armZZos[24]= - armZZos[12]*armZZos[21];
   armZZos[26]= - armZZos[11]*armZZos[31];
   armZZos[24]=armZZos[24] + armZZos[26];
   armZZos[26]=armZZos[16]*armZZos[30];
   armZZos[24]=4./3.*armZZos[24] + armZZos[26];
   armZZos[24]=mmZ*armZZos[24];
   armZZos[23]=armZZos[23]*mmt;
   armZZos[26]= - armZZos[23] + armZZos[31];
   armZZos[26]=armZZos[13]*mmt*armZZos[26];
   armZZos[23]= - 4./3.*armZZos[23] + armZZos[29];
   armZZos[23]=armZZos[8]*armZZos[23];
   armZZos[21]=armZZos[10]*armZZos[21];
   armZZos[19]=armZZos[22] + armZZos[23] + 16./3.*armZZos[26] + 
   armZZos[24] - 2*armZZos[21] + 1./3.*armZZos[20] + armZZos[19] + 4*
   armZZos[25];
   armZZos[19]=armZZos[1]*armZZos[19];
   armZZos[17]=11./9.*armZZos[18] + armZZos[17] - 20./9.;
   armZZos[18]=mmZ*armZZos[12];
   armZZos[18]= - 8./3.*armZZos[18] - 31./3. - 4*armZZos[10];
   armZZos[17]=armZZos[15]*armZZos[17]*armZZos[18];

      mZZosret = armZZos[17] + armZZos[19];
      return mZZosret;
}
