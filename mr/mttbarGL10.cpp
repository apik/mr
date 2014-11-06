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
tt::xgl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[11], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbarGL[3]=pow(SW,-1);
    armttbarGL[4]=pow(MMW,-1);
    armttbarGL[5]=Tsil::A(MMH,mu2);
    armttbarGL[6]=Tsil::A(MMt,mu2);
    armttbarGL[7]=pow(MMH,-1);
    armttbarGL[8]=std::real(Tsil::B(0,0,MMt,mu2));
   armttbarGL[9]=1./2.*armttbarGL[2];
   armttbarGL[10]=armttbarGL[6]*armttbarGL[7];
   armttbarGL[10]=1./8.*armttbarGL[8] + armttbarGL[9] - 3*
   armttbarGL[10];
   armttbarGL[10]=MMt*armttbarGL[10];
   armttbarGL[9]= - MMH*armttbarGL[9];
   armttbarGL[9]=armttbarGL[9] + armttbarGL[5] + armttbarGL[6];
   armttbarGL[9]=1./4.*armttbarGL[9] + armttbarGL[10];

      mttbarGLret = armttbarGL[9]*armttbarGL[4]*pow(armttbarGL[3],2)*
      armttbarGL[1];
      return mttbarGLret;
}
