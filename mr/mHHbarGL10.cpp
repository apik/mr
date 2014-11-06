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
HH<OS>::xgl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbarGL[12], mHHbarGLret;

    armHHbarGL[1]=double(boson);
    armHHbarGL[2]=Tsil::B(MMH,MMH,MMH,mu2);
    armHHbarGL[3]=pow(SW,-1);
    armHHbarGL[4]=pow(MMW,-1);
    armHHbarGL[5]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbarGL[6]=pow(MMH,-1);
    armHHbarGL[7]=Tsil::A(MMH,mu2);
    armHHbarGL[8]=Tsil::A(MMt,mu2);
    armHHbarGL[9]=std::real(Tsil::B(0,0,MMH,mu2));
   armHHbarGL[10]= - MMt*armHHbarGL[5];
   armHHbarGL[10]=armHHbarGL[10] - armHHbarGL[8];
   armHHbarGL[11]=2*armHHbarGL[6];
   armHHbarGL[10]=armHHbarGL[11]*armHHbarGL[10];
   armHHbarGL[10]=1./2.*armHHbarGL[5] + armHHbarGL[10];
   armHHbarGL[10]=MMt*armHHbarGL[10];
   armHHbarGL[11]=3./8.*armHHbarGL[2] + 1./8.*armHHbarGL[9];
   armHHbarGL[11]=MMH*armHHbarGL[11];
   armHHbarGL[10]=1./4.*armHHbarGL[7] + armHHbarGL[10] + armHHbarGL[11]
   ;

      mHHbarGLret = 3*armHHbarGL[10]*armHHbarGL[4]*pow(armHHbarGL[3],2)
      *armHHbarGL[1];
      return mHHbarGLret;
}
