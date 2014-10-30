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
ZZ<OS>::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZbar[32], mZZbarret;

    armZZbar[1]=double(nH);
    armZZbar[2]=pow(CW,-1);
    armZZbar[3]=pow(MMH,-1);
    armZZbar[4]=pow(MMZ,-1);
    armZZbar[5]=pow(SW,-1);
    armZZbar[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    armZZbar[7]=Tsil::A(MMt,mu2);
    armZZbar[8]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    armZZbar[9]=Tsil::Aeps(MMt,mu2);
    armZZbar[10]=std::real(Tsil::B(0,0,MMZ,mu2));
    armZZbar[11]=prottttt0->M(0);
    armZZbar[12]=prot00000->M(0);
    armZZbar[13]=prottttt0->Vzxyv(0);
    armZZbar[14]=prottttt0->Suxv(0);
    armZZbar[15]=double(nL);
    armZZbar[16]=1/(4*MMt - MMZ);
   armZZbar[17]=pow(armZZbar[5],2);
   armZZbar[18]=pow(armZZbar[2],2);
   armZZbar[19]= - 7./9.*armZZbar[18] + armZZbar[17] + 64./9.;
   armZZbar[20]=armZZbar[19]*armZZbar[6];
   armZZbar[21]=29./3.*armZZbar[18] - 128./3. - armZZbar[17];
   armZZbar[21]=2*armZZbar[21] - armZZbar[20];
   armZZbar[21]=armZZbar[6]*armZZbar[21];
   armZZbar[21]=armZZbar[21] - 299./9.*armZZbar[18] - 64./9. - 35*
   armZZbar[17];
   armZZbar[22]=2./3.*armZZbar[13] + 1./3.*armZZbar[11];
   armZZbar[22]=armZZbar[19]*armZZbar[22];
   armZZbar[23]=armZZbar[17] + armZZbar[18];
   armZZbar[24]=armZZbar[23]*armZZbar[3];
   armZZbar[22]=8*armZZbar[24] + armZZbar[22];
   armZZbar[25]=4*MMt;
   armZZbar[22]=armZZbar[22]*armZZbar[25];
   armZZbar[21]=1./3.*armZZbar[21] + armZZbar[22];
   armZZbar[21]=MMt*armZZbar[21];
   armZZbar[22]=7*armZZbar[17];
   armZZbar[26]=95./9.*armZZbar[18] - 128./9. + armZZbar[22];
   armZZbar[26]=armZZbar[9]*armZZbar[26];
   armZZbar[27]=17./9.*armZZbar[18] + armZZbar[17] - 32./9.;
   armZZbar[28]= - armZZbar[14]*armZZbar[27];
   armZZbar[26]=1./3.*armZZbar[26] + armZZbar[28];
   armZZbar[28]=armZZbar[6] + 2;
   armZZbar[28]=armZZbar[28]*armZZbar[23];
   armZZbar[29]=armZZbar[6]*armZZbar[28];
   armZZbar[29]=armZZbar[29] - armZZbar[23];
   armZZbar[29]=MMt*armZZbar[29];
   armZZbar[30]=2*armZZbar[9] - armZZbar[14];
   armZZbar[30]=armZZbar[23]*armZZbar[30];
   armZZbar[29]=armZZbar[29] + armZZbar[30];
   armZZbar[30]=armZZbar[4]*MMt;
   armZZbar[29]=armZZbar[29]*armZZbar[30];
   armZZbar[21]=2./3.*armZZbar[29] + 2*armZZbar[26] + armZZbar[21];
   armZZbar[26]=2*armZZbar[4];
   armZZbar[21]=armZZbar[21]*armZZbar[26];
   armZZbar[26]=25./9.*armZZbar[18] + armZZbar[17] - 64./9.;
   armZZbar[29]=armZZbar[26]*armZZbar[6];
   armZZbar[22]= - armZZbar[29] + 143./9.*armZZbar[18] - 320./9. + 
   armZZbar[22];
   armZZbar[22]=armZZbar[6]*armZZbar[22];
   armZZbar[23]= - armZZbar[11]*armZZbar[23];
   armZZbar[31]=armZZbar[13]*armZZbar[27];
   armZZbar[23]=armZZbar[23] - 4./3.*armZZbar[31];
   armZZbar[23]=armZZbar[23]*armZZbar[25];
   armZZbar[25]= - armZZbar[29] + 25./3.*armZZbar[18] - 64./3. + 3*
   armZZbar[17];
   armZZbar[25]=MMZ*armZZbar[6]*armZZbar[25];
   armZZbar[29]=armZZbar[9]*armZZbar[26];
   armZZbar[25]= - 4*armZZbar[29] + armZZbar[25];
   armZZbar[25]=armZZbar[16]*armZZbar[25];
   armZZbar[29]= - armZZbar[16]*MMZ;
   armZZbar[29]=armZZbar[29] - 1;
   armZZbar[29]=armZZbar[26]*armZZbar[29];
   armZZbar[19]=armZZbar[19]*armZZbar[30];
   armZZbar[19]=4./3.*armZZbar[19] + armZZbar[29];
   armZZbar[19]=armZZbar[8]*armZZbar[19];
   armZZbar[28]=armZZbar[28]*armZZbar[30];
   armZZbar[20]=armZZbar[20] - armZZbar[28];
   armZZbar[20]=211./9.*armZZbar[18] - 448./9. + 11*armZZbar[17] - 2*
   armZZbar[20];
   armZZbar[24]=armZZbar[7]*armZZbar[24];
   armZZbar[20]= - 12*armZZbar[24] + 1./3.*armZZbar[20];
   armZZbar[20]=armZZbar[4]*armZZbar[20];
   armZZbar[24]= - armZZbar[6] + 1;
   armZZbar[24]=armZZbar[16]*armZZbar[26]*armZZbar[24];
   armZZbar[20]=armZZbar[24] + armZZbar[20];
   armZZbar[20]=armZZbar[7]*armZZbar[20];
   armZZbar[24]=5./9.*armZZbar[18] + armZZbar[17] - 8./9.;
   armZZbar[26]=armZZbar[12]*MMZ;
   armZZbar[28]=2*armZZbar[10] + 4./3.*armZZbar[26];
   armZZbar[24]=armZZbar[24]*armZZbar[28];
   armZZbar[28]=937./18.*armZZbar[18] - 860./9. + 77./2.*armZZbar[17];
   armZZbar[29]=4./3.*MMZ;
   armZZbar[27]=armZZbar[11]*armZZbar[27]*armZZbar[29];
   armZZbar[19]=4*armZZbar[20] + armZZbar[19] + armZZbar[25] + 
   armZZbar[21] + armZZbar[23] + armZZbar[27] + 1./3.*armZZbar[28] + 
   armZZbar[22] + armZZbar[24];
   armZZbar[19]=armZZbar[1]*armZZbar[19];
   armZZbar[17]=11./9.*armZZbar[18] + armZZbar[17] - 20./9.;
   armZZbar[18]=31 + 8*armZZbar[26];
   armZZbar[18]=1./3.*armZZbar[18] + 4*armZZbar[10];
   armZZbar[17]=armZZbar[15]*armZZbar[17]*armZZbar[18];

      mZZbarret = armZZbar[17] + armZZbar[19];
      return mZZbarret;
}
