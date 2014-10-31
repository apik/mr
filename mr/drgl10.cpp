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
dr::drgl10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> ardrgl[9], drglret;

    ardrgl[1]=double(boson);
    ardrgl[2]=pow(SW,-1);
    ardrgl[3]=pow(MMW,-1);
    ardrgl[4]=Tsil::A(MMH,mu2);
    ardrgl[5]=Tsil::A(MMt,mu2);
    ardrgl[6]=pow(MMH,-1);
   ardrgl[7]=ardrgl[4] + MMt;
   ardrgl[8]=ardrgl[6]*MMt;
   ardrgl[8]=1./2. - 2*ardrgl[8];
   ardrgl[8]=ardrgl[5]*ardrgl[8];
   ardrgl[7]=ardrgl[8] + 1./4.*ardrgl[7];
   ardrgl[7]=3*ardrgl[7] - 1./8.*MMH;

      drglret = ardrgl[7]*ardrgl[3]*pow(ardrgl[2],2)*ardrgl[1];
      return drglret;
}
