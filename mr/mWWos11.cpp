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

#include <WW.hpp>
std::complex<long double>
WW<MS>::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armWWos[23], mWWosret;

    armWWos[1]=double(nH);
    armWWos[2]=pow(mmZ,-1);
    armWWos[3]=pow(mmH,-1);
    armWWos[4]=pow(s,-1);
    armWWos[5]=pow(c,-1);
    armWWos[6]=Tsil::B(0,mmt,mmW,mu2);
    armWWos[7]=Tsil::A(mmt,mu2);
    armWWos[8]=pow(mmt,-1);
    armWWos[9]=Tsil::Aeps(mmt,mu2);
    armWWos[10]=prot00tt0->M(0);
    armWWos[11]=prot00tt0->Tuxv(0);
    armWWos[12]=double(nL);
    armWWos[13]=std::real(Tsil::B(0,0,mmW,mu2));
    armWWos[14]=prot00000->M(0);
   armWWos[15]= - 7 + 5./3.*armWWos[6];
   armWWos[16]=2*armWWos[6];
   armWWos[15]=armWWos[15]*armWWos[16];
   armWWos[16]=mmt*armWWos[3];
   armWWos[17]=armWWos[7]*armWWos[3];
   armWWos[15]=8./3.*armWWos[11] + armWWos[15] + 16*armWWos[17] - 17./6.
    + 48*armWWos[16];
   armWWos[15]=armWWos[15]*mmt;
   armWWos[16]= - 23 + 14*armWWos[6];
   armWWos[16]=12*armWWos[17] + 1./3.*armWWos[16];
   armWWos[17]=4*armWWos[7];
   armWWos[16]=armWWos[16]*armWWos[17];
   armWWos[17]=pow(armWWos[7],2);
   armWWos[18]=armWWos[17]*armWWos[8];
   armWWos[15]=8*armWWos[18] + 8./3.*armWWos[9] + armWWos[15] + 
   armWWos[16];
   armWWos[16]=pow(armWWos[5],2);
   armWWos[18]=armWWos[16]*armWWos[1];
   armWWos[19]= - armWWos[15]*armWWos[18];
   armWWos[20]= - 11./3. + 6*armWWos[6];
   armWWos[20]=armWWos[20]*armWWos[7];
   armWWos[21]=mmt*armWWos[10];
   armWWos[22]= - 2./3.*armWWos[21] + 5*armWWos[6];
   armWWos[22]=armWWos[22]*mmt;
   armWWos[20]= - 4./3.*armWWos[9] + armWWos[20] - armWWos[22];
   armWWos[20]=mmt*armWWos[20];
   armWWos[22]=4./3.*armWWos[11];
   armWWos[22]=armWWos[22]*pow(mmt,2);
   armWWos[17]= - armWWos[22] + 20./3.*armWWos[17] + armWWos[20];
   armWWos[20]=2*armWWos[2];
   armWWos[17]=armWWos[17]*armWWos[20];
   armWWos[18]= - armWWos[1] - armWWos[18];
   armWWos[16]=armWWos[16]*armWWos[18]*armWWos[17];
   armWWos[16]=armWWos[19] + armWWos[16];
   armWWos[16]=armWWos[2]*armWWos[16];
   armWWos[15]= - armWWos[15] - armWWos[17];
   armWWos[15]=armWWos[2]*armWWos[15];
   armWWos[17]= - 1 - 2./3.*armWWos[6];
   armWWos[17]=armWWos[6]*armWWos[17];
   armWWos[17]=armWWos[21] + armWWos[17];
   armWWos[18]=1 - armWWos[6];
   armWWos[18]=armWWos[7]*armWWos[18];
   armWWos[18]=armWWos[18] - armWWos[9];
   armWWos[18]=16./3.*armWWos[18];
   armWWos[18]=armWWos[8]*armWWos[18];
   armWWos[15]=armWWos[15] - 16./3.*armWWos[11] - 13 + armWWos[18] + 4*
   armWWos[17];
   armWWos[15]=armWWos[1]*armWWos[15];
   armWWos[17]=armWWos[14]*armWWos[12];
   armWWos[18]=armWWos[1]*armWWos[10];
   armWWos[17]=armWWos[17] + armWWos[18];
   armWWos[18]=8./3.*mmZ;
   armWWos[17]=armWWos[17]*armWWos[18];
   armWWos[18]= - 4*armWWos[13] - 31./3.;
   armWWos[18]=armWWos[12]*armWWos[18];
   armWWos[15]= - armWWos[17] + armWWos[18] + armWWos[15];
   armWWos[15]=armWWos[15]*pow(armWWos[4],2);

      mWWosret = armWWos[15] + armWWos[16] + armWWos[17];
      return mWWosret;
}
