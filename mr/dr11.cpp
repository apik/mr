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

#include <dr.hpp>
std::complex<long double>
dr::dr11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardr[13], drret;

    ardr[1]=double(nH);
    ardr[2]=pow(CW,-1);
    ardr[3]=pow(MMH,-1);
    ardr[4]=pow(MMZ,-1);
    ardr[5]=pow(SW,-1);
    ardr[6]=Tsil::I2(0,0,MMt,mu2);
    ardr[7]=Tsil::A(MMt,mu2);
    ardr[8]=pow(MMt,-1);
    ardr[9]=Tsil::Aeps(MMt,mu2);
   ardr[10]= - ardr[8] + 12*ardr[3];
   ardr[11]=2*ardr[7];
   ardr[10]=ardr[10]*ardr[11];
   ardr[10]=ardr[10] - 5;
   ardr[10]=ardr[10]*ardr[11];
   ardr[11]=ardr[9] - ardr[6];
   ardr[12]=64*ardr[3];
   ardr[12]=ardr[12]*pow(MMt,2);
   ardr[10]=ardr[10] - ardr[12] + 37./2.*MMt - 4*ardr[11];
   ardr[11]= - pow(ardr[2],2);
   ardr[12]= - pow(ardr[5],2);
   ardr[11]=ardr[11] + ardr[12];

      drret = ardr[11]*ardr[10]*ardr[4]*ardr[1];
      return drret;
}
