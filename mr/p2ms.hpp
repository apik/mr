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

#ifndef __P2MS_HPP__
#define __P2MS_HPP__

#include "sminput.hpp"

#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp"
#include "bb.hpp"
#include "dr.hpp"


class P2MS
{
  OSinput       oi;
    
  long double  aEW;
  long double aQCD;
  long double   Gf;
  long double   mu;
  

  bb*      bp;
  WW<OS>*  wp;
  ZZ<OS>*  zp;
  HH<OS>*  hp;
  tt<OS>*  tp;
  dr<OS>* drp;
  
  
public:
  P2MS(const OSinput & oi_,
       const long double &  Gf_ = pdg2014::Gf,
       const long double &  as_ = pdg2014::asMZ,
       const long double &  mu_ = pdg2014::MZ );

  long double   a1(size_t nL = 2, size_t nH = 1);
  long double   a2(size_t nL = 2, size_t nH = 1);
  long double   as(size_t nL = 2, size_t nH = 1);
  long double   at(size_t nL = 2, size_t nH = 1);
  long double   ab(size_t nL = 2, size_t nH = 1);
  long double alam(size_t nL = 2, size_t nH = 1);

  
  long double  g1(size_t nL = 2, size_t nH = 1);
  long double  g2(size_t nL = 2, size_t nH = 1);
  long double  gs(size_t nL = 2, size_t nH = 1);
  long double  yt(size_t nL = 2, size_t nH = 1);
  long double  yb(size_t nL = 2, size_t nH = 1);
  long double lam(size_t nL = 2, size_t nH = 1);
  long double mu0(size_t nL = 2, size_t nH = 1);
  long double vev(size_t nL = 2, size_t nH = 1);
  
  MSinput getMSpar();
};


#endif  // __P2MS_HPP__
