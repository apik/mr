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

#include <bb.hpp>
std::complex<long double>
bb::mgl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbarGL[9], mbbbarGLret;

    armbbbarGL[1]=double(boson);
    armbbbarGL[2]=pow(SW,-1);
    armbbbarGL[3]=pow(MMW,-1);
    armbbbarGL[4]=Tsil::A(MMH,mu2);
    armbbbarGL[5]=Tsil::A(MMt,mu2);
    armbbbarGL[6]=pow(MMH,-1);
   armbbbarGL[7]=armbbbarGL[4] + armbbbarGL[5];
   armbbbarGL[8]=armbbbarGL[6]*armbbbarGL[5];
   armbbbarGL[8]=1./16. - 3*armbbbarGL[8];
   armbbbarGL[8]=MMt*armbbbarGL[8];
   armbbbarGL[7]=3./8.*armbbbarGL[7] + armbbbarGL[8];

      mbbbarGLret = armbbbarGL[7]*armbbbarGL[3]*pow(armbbbarGL[2],2)*
      armbbbarGL[1];
      return mbbbarGLret;
}
