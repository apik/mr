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


class P2MSnLnH
{
  OSinput       oi;
    
  long double  aEW;
  long double aQCD;
  long double   Gf;
  long double   mu;
  

  bb<OS>*  bp;
  WW<OS>*  wp;
  ZZ<OS>*  zp;
  HH<OS>*  hp;
  tt<OS>*  tp;
  dr<OS>* drp;
  
  
public:
  P2MSnLnH(const OSinput & oi_,
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

  std::vector<long double> runningCouplings();
};


// SM version with fixed nL=2, nH=1
class P2MS
{
  OSinput       oi;
    
  long double  aEW;
  long double aQCD;
  long double   Gf;
  long double   mu;
  

  bb<OS>*  bp;
  WW<OS>*  wp;
  ZZ<OS>*  zp;
  HH<OS>*  hp;
  tt<OS>*  tp;
  dr<OS>* drp;

  // corrections
  long double dbplus1;
  long double dWplus1;
  long double dZplus1;
  long double dHplus1;
  long double dtplus1;
  long double dRplus1;

  unsigned ord;

public:
  P2MS(const OSinput & oi_,
       const long double &  Gf_ = pdg2014::Gf,
       const long double &  as_ = pdg2014::asMZ,
       const long double &  mu_ = pdg2014::MZ,
       unsigned ord_ = order::x01|order::x10|order::x02|order::x11|order::x20|order::x03 );
  
  long double   a1() const;
  long double   a2() const;
  long double   as() const;
  long double   at() const;
  long double   ab() const;
  long double alam() const;
  
  
  long double  g1() const;
  long double  g2() const;
  long double  gs() const;
  long double  yt() const;
  long double  yb() const;
  long double lam() const;
  long double mu0() const;
  long double vev() const;
  
  MSinput getMSpar();

  std::vector<long double> runningCouplings() const;

  std::vector<long double> ai() const;
    
  long double scale() const
  {
    return mu;
  };
  
};


#endif  // __P2MS_HPP__
