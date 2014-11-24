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

#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__
const long double    Pi = M_PI;
const long double Zeta2 = (Pi*Pi)/6.;
const long double Zeta3 = 1.20205690315959428539973816151144999076498629234049888179227155534L;
const long double Zeta4 = (Pi*Pi*Pi*Pi)/90.;
const long double Zeta5 = 1.03692775514336992633136548645703416805708091950191281197419267790L;

// a4 = Li[4,1/2]
const long double a4    = 0.5174790616738993863307581618988629456223774751413792582443193479770L;
// a5 = Li[5,1/2]
const long double a5    = 0.5084005792422687074591088492585899413195411256648216487244977963526L;

const long double S1    = Pi/sqrt(3.);
const long double S2    = 0.26043413763216209895572914320803078545504477884842847340736668765L;

// g5 = -i/24*EPAIR*eps(mu1,mu2,mu3,mu4)*g(mu1)*g(mu2)*g(mu3)*g(mu4)
// Sp[g5*g5] = 4
const long double EPAIR2 = -1.;


// Latest input from:
// Review of Particle Physics, 2014
// K.A. Olive et al. (Particle Data Group), Chin. Phys. C, 38, 090001 (2014).
namespace pdg2014
{
  const double   MW = 80.385;
  const double   MZ = 91.1876;
  const double   MH = 125.7;
  const double   Mt = 173.21;

  // MS running
  const double   mb = 4.18;
  // Mb(mb(mu=mb)) in QCD 
  const double   Mb = 4.93482;

  const double   Gf = 1.1663787e-5;
  
  const double  aMZ = 1./127.940;
  const double asMZ = 0.1185;
}


// Review of Particle Physics, 2012
// J. Beringer et al. (Particle Data Group), Phys. Rev. D86, 010001 (2012).
namespace pdg2012
{
  const double   MW = 80.385;
  const double   MZ = 91.1876;
  // double   MH = ;
  const double   Mt = 173.5;

  // MS running
  const double   mb = 4.18;
  // Mb(mb(mu=mb)) in QCD 
  const double   Mb = 4.93482;
  
  const double   Gf = 1.1663787e-5;
  
  const double  aMZ = 1./127.944;
  const double asMZ = 0.1184;
}


// Review of Particle Physics, 2010
// K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010).
namespace pdg2010
{
  const double   MW = 80.399;
  const double   MZ = 91.1876;
  // double   MH = ;          
  const double   Mt = 172.0;

  // MS running
  const double   mb = 4.19;
  // Mb(mb(mu=mb)) in QCD using 3-loop 
  const double   Mb = 4.9456;

  const double   Gf = 1.16637e-5;
  
  const double  aMZ = 1./127.916;
  const double asMZ = 0.1184;
}


#endif  // __CONSTANTS_HPP__
