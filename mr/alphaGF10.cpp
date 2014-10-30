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

#include <alphaGF.hpp>
std::complex<long double>
alphaGF::a10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aralphaGF[35], alphaGFret;

    aralphaGF[1]=double(nL + nH);
    aralphaGF[2]=std::real(Tsil::B(0,0,MMZ,mu2));
    aralphaGF[3]=pow(SW,-1);
    aralphaGF[4]=std::real(Tsil::B(0,0,MMW,mu2));
    aralphaGF[5]=double(nH);
    aralphaGF[6]=pow(CW,-1);
    aralphaGF[7]=pow(MMZ,-1);
    aralphaGF[8]=Tsil::B(MMt,MMt,MMZ,mu2);
    aralphaGF[9]=Tsil::B(0,MMt,MMW,mu2);
    aralphaGF[10]=Tsil::A(MMt,mu2);
    aralphaGF[11]=double(nL);
    aralphaGF[12]=double(boson);
    aralphaGF[13]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aralphaGF[14]=Tsil::B(MMW,MMH,MMW,mu2);
    aralphaGF[15]=Tsil::B(MMW,MMZ,MMW,mu2);
    aralphaGF[16]=Tsil::B(MMW,MMW,MMZ,mu2);
    aralphaGF[17]=Tsil::A(MMH,mu2);
    aralphaGF[18]=Tsil::A(MMZ,mu2);
    aralphaGF[19]=Tsil::A(MMW,mu2);
    aralphaGF[20]=1/( - MMW + MMH);
   aralphaGF[21]=pow(aralphaGF[6],2);
   aralphaGF[22]=1./4.*aralphaGF[21];
   aralphaGF[23]=17 + aralphaGF[21];
   aralphaGF[23]=aralphaGF[23]*aralphaGF[22];
   aralphaGF[24]=pow(aralphaGF[3],2);
   aralphaGF[25]=2*aralphaGF[24];
   aralphaGF[26]= - 31./4. + aralphaGF[25];
   aralphaGF[26]=aralphaGF[24]*aralphaGF[26];
   aralphaGF[23]=aralphaGF[23] + aralphaGF[26];
   aralphaGF[23]=aralphaGF[18]*aralphaGF[23];
   aralphaGF[26]= - 1./8. - aralphaGF[14];
   aralphaGF[26]=aralphaGF[21]*aralphaGF[26];
   aralphaGF[27]=aralphaGF[24] - 1;
   aralphaGF[28]=aralphaGF[14]*aralphaGF[27];
   aralphaGF[28]= - 1./8. + aralphaGF[28];
   aralphaGF[28]=aralphaGF[28]*aralphaGF[24];
   aralphaGF[29]=pow(aralphaGF[3],4);
   aralphaGF[30]=aralphaGF[29]*aralphaGF[13];
   aralphaGF[26]= - aralphaGF[30] + aralphaGF[28] + aralphaGF[26];
   aralphaGF[26]=MMH*aralphaGF[26];
   aralphaGF[23]=aralphaGF[26] + aralphaGF[23];
   aralphaGF[26]=aralphaGF[13] - aralphaGF[14];
   aralphaGF[26]=aralphaGF[29]*aralphaGF[26];
   aralphaGF[28]=pow(aralphaGF[6],4);
   aralphaGF[31]=aralphaGF[14]*aralphaGF[28];
   aralphaGF[26]=aralphaGF[31] + aralphaGF[26];
   aralphaGF[26]=MMH*aralphaGF[26];
   aralphaGF[31]=aralphaGF[28] - aralphaGF[29];
   aralphaGF[32]= - aralphaGF[19]*aralphaGF[31];
   aralphaGF[29]= - aralphaGF[18]*aralphaGF[29];
   aralphaGF[28]=aralphaGF[17]*aralphaGF[28];
   aralphaGF[26]=aralphaGF[32] + aralphaGF[28] + aralphaGF[29] + 
   aralphaGF[26];
   aralphaGF[26]=aralphaGF[7]*MMH*aralphaGF[26];
   aralphaGF[28]=1./2.*aralphaGF[21];
   aralphaGF[29]=9 - 1./6.*aralphaGF[21];
   aralphaGF[29]=aralphaGF[29]*aralphaGF[28];
   aralphaGF[25]=19./2. - aralphaGF[25];
   aralphaGF[25]=aralphaGF[25]*aralphaGF[24];
   aralphaGF[25]=1./3.*aralphaGF[25] - 4 + aralphaGF[29];
   aralphaGF[25]=aralphaGF[19]*aralphaGF[25];
   aralphaGF[29]= - aralphaGF[21] - aralphaGF[24];
   aralphaGF[29]=aralphaGF[17]*aralphaGF[29];
   aralphaGF[23]=1./12.*aralphaGF[26] + aralphaGF[25] + 1./3.*
   aralphaGF[23] + 1./4.*aralphaGF[29];
   aralphaGF[23]=aralphaGF[7]*aralphaGF[23];
   aralphaGF[25]= - 2 + 3./4.*aralphaGF[24];
   aralphaGF[26]=11*aralphaGF[24];
   aralphaGF[25]=aralphaGF[25]*aralphaGF[26];
   aralphaGF[26]=4 + aralphaGF[22];
   aralphaGF[26]=aralphaGF[26]*aralphaGF[21];
   aralphaGF[26]=aralphaGF[25] + 8 + 1./3.*aralphaGF[26];
   aralphaGF[26]=aralphaGF[15]*aralphaGF[26];
   aralphaGF[29]=aralphaGF[24] - 2;
   aralphaGF[32]= - aralphaGF[14]*aralphaGF[29];
   aralphaGF[33]=aralphaGF[17] - aralphaGF[19];
   aralphaGF[33]=aralphaGF[20]*aralphaGF[33];
   aralphaGF[32]=3./4.*aralphaGF[33] + 63./8. + aralphaGF[32];
   aralphaGF[32]=aralphaGF[24]*aralphaGF[32];
   aralphaGF[33]=pow(CW,2);
   aralphaGF[33]=4*aralphaGF[33];
   aralphaGF[25]= - aralphaGF[33] - 41./3. - aralphaGF[25];
   aralphaGF[25]=aralphaGF[16]*aralphaGF[25];
   aralphaGF[34]= - 44 - 1./8.*aralphaGF[21];
   aralphaGF[23]=aralphaGF[23] + aralphaGF[26] + aralphaGF[25] + 
   aralphaGF[30] - aralphaGF[33] + 1./3.*aralphaGF[34] + aralphaGF[32];
   aralphaGF[23]=aralphaGF[12]*aralphaGF[23];
   aralphaGF[25]=aralphaGF[27]*aralphaGF[24];
   aralphaGF[21]= - aralphaGF[21] + aralphaGF[25];
   aralphaGF[21]=aralphaGF[9]*aralphaGF[21];
   aralphaGF[25]= - MMt*aralphaGF[9];
   aralphaGF[25]= - aralphaGF[10] + aralphaGF[25];
   aralphaGF[25]=aralphaGF[7]*aralphaGF[31]*aralphaGF[25];
   aralphaGF[21]=aralphaGF[21] + aralphaGF[25];
   aralphaGF[25]=1./2.*aralphaGF[24];
   aralphaGF[26]= - 8./3. - aralphaGF[25];
   aralphaGF[26]=aralphaGF[26]*aralphaGF[24];
   aralphaGF[26]=32./9. + aralphaGF[26];
   aralphaGF[26]=aralphaGF[8]*aralphaGF[26];
   aralphaGF[21]=aralphaGF[26] - 29./12.*aralphaGF[24] + 32./9. + 
   aralphaGF[22] + 1./2.*aralphaGF[21];
   aralphaGF[21]=MMt*aralphaGF[21];
   aralphaGF[22]= - 19./6.*aralphaGF[24] + 32./9. - aralphaGF[28];
   aralphaGF[22]=aralphaGF[10]*aralphaGF[22];
   aralphaGF[21]=aralphaGF[22] + aralphaGF[21];
   aralphaGF[21]=aralphaGF[7]*aralphaGF[21];
   aralphaGF[22]= - 4./3. + aralphaGF[25];
   aralphaGF[22]=aralphaGF[22]*aralphaGF[24];
   aralphaGF[22]=16./9. + aralphaGF[22];
   aralphaGF[22]=aralphaGF[8]*aralphaGF[22];
   aralphaGF[26]=aralphaGF[29]*aralphaGF[24];
   aralphaGF[27]= - aralphaGF[9]*aralphaGF[26];
   aralphaGF[25]= - 2./3. + aralphaGF[25];
   aralphaGF[24]=aralphaGF[25]*aralphaGF[24];
   aralphaGF[24]=4./9. + aralphaGF[24];
   aralphaGF[24]=aralphaGF[2]*aralphaGF[24];
   aralphaGF[21]=aralphaGF[21] + aralphaGF[24] + aralphaGF[27] - 20./27.
    + aralphaGF[22];
   aralphaGF[21]=aralphaGF[5]*aralphaGF[21];
   aralphaGF[22]=20./9. + aralphaGF[26];
   aralphaGF[22]=aralphaGF[11]*aralphaGF[22];
   aralphaGF[24]=4 + aralphaGF[26];
   aralphaGF[25]=1./3.*aralphaGF[1];
   aralphaGF[24]=aralphaGF[24]*aralphaGF[25];
   aralphaGF[22]=aralphaGF[22] + aralphaGF[24];
   aralphaGF[22]=aralphaGF[2]*aralphaGF[22];
   aralphaGF[24]=aralphaGF[4]*aralphaGF[26];
   aralphaGF[26]= - 20./27. - aralphaGF[24];
   aralphaGF[26]=aralphaGF[11]*aralphaGF[26];
   aralphaGF[24]= - 4./3. - aralphaGF[24];
   aralphaGF[24]=aralphaGF[24]*aralphaGF[25];

      alphaGFret = aralphaGF[21] + aralphaGF[22] + aralphaGF[23] + 
      aralphaGF[24] + aralphaGF[26];
      return alphaGFret;
}
