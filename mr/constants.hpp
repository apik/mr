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
const long double Zeta3 = 1.20205690315959428539973816151144999076498629234049888179227155534;
const long double Zeta4 = (Pi*Pi*Pi*Pi)/90.;
const long double Zeta5 = 1.03692775514336992633136548645703416805708091950191281197419267790;

const long double S1    = Pi/sqrt(3.);
const long double S2    = 0.26043413763216209895572914320803078545504477884842847340736668765;

// g5 = -i/24*EPAIR*eps(mu1,mu2,mu3,mu4)*g(mu1)*g(mu2)*g(mu3)*g(mu4)
// Sp[g5*g5] = 4
const long double EPAIR2 = -1.;

#endif  // __CONSTANTS_HPP__
