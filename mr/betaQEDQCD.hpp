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

#ifndef __BETAQEDQCD_HPP_
#define __BETAQEDQCD_HPP_

// #include <boost/numeric/odeint/integrator_adaptive_stepsize.hpp>
#include <stdexcept>

#include "sminput.hpp"
#include "constants.hpp"
#include "betaQCD.hpp"


typedef std::vector<long double> state_type;

class BetaQEDQCD
{
  
  BetaQCD* pbQCD;
  
  size_t ng;
  
  bool MultiplyByMinus1;
public:

  BetaQEDQCD( size_t ng_, bool MultiplyByMinus1_ = false) : ng(ng_), MultiplyByMinus1(MultiplyByMinus1_)
  {
    pbQCD = new BetaQCD(4, 2*ng, MultiplyByMinus1);
  }
  
  void operator() (const state_type &av4pi, state_type &dadt, const double t)
  {
    
    long double NG = static_cast<long double>(ng);
    
    double minusC = MultiplyByMinus1 ? -1. : 1.;

    if(av4pi.size() != 2)
      throw std::logic_error("ERROR: for mu^2 running 2 constants in input needed: aQCD and aQED");

    pbQCD->operator()(av4pi, dadt, t);

    long double as4pi = av4pi[0], a4pi = av4pi[1];

    // 1-loop alphaEM beta-function in full SM with
    // QCD corrections upto five loops from:
    // Baikov, Chetyrkin, Kuhn, Rittinger, JHEP 1207 (2012) 017
    // arXiv:1206.1284 [hep-ph]
    dadt[1] =                   // We need -\\beta if evolve from t=0 to t < 0
      minusC*(
              pow(a4pi,2) *     // 1-loop pure EW
              ( -7. + 32/9. * NG + (80/9.) *
                as4pi * NG +    // 2-loop QCD
                pow(as4pi,2) *  // 3-loop QCD
                ((2500/27.) * (NG) + (-440/27.) * (pow(NG,2))) + 
                pow(as4pi,3) *  // 4-loop QCD
                ((209740/243.) * (NG) + (-83080/243.) * (pow(NG,2)) 
                 + (-6160/243.) * (pow(NG,3)) + (35200/81.) * ((NG) * (Zeta3)) + 
                 (-12160/27.) * ((pow(NG,2)) * (Zeta3))) +
                pow(as4pi,4) *  // 5-loop QCD
                ((13326745/1458.) * (NG) + (-1797080/243.) * (pow(NG,2)) + 
                 (-60490/243.) * (pow(NG,3)) + (2140/81.) * (pow(NG,4)) + 
                 (7293400/243.) * ((NG) * (Zeta3)) + (-7313200/243.) * ((pow(NG,2)) * (Zeta3)) +
                 (262880/81.) * ((pow(NG,3)) * (Zeta3)) + (320/9.) * 
                 ((pow(NG,4)) * (Zeta3)) + (-48400/9.) * ((NG) * (Zeta4)) +
                 (58960/9.)  * ((pow(NG,2)) * (Zeta4)) +
                 (-3040/3.) * ((pow(NG,3)) * (Zeta4)) + (-1255000/81.) * ((NG) * (Zeta5)) +
                 (3322000/243.) * ((pow(NG,2)) * (Zeta5)) + (-1600/27.) * ((pow(NG,3)) * (Zeta5))))
              );
    
  }
};



std::pair<double,double> runQEDQCD(long double aStart, long double asStart, long double muStart, long double muEnd, size_t ng = 3);



class AlphaQEDQCD
{
  double muStart;
  double  aStart;
  double asStart;
  size_t   loops;
  size_t      ng;
  OSinput     oi;

public:

  // Running with fixed nf
  AlphaQEDQCD(long double aMZ = pdg2014::aMZ,long double asMZ = pdg2014::asMZ, long double mu = pdg2014::MZ, size_t ng_ = 3) : aStart(aMZ), asStart(asMZ), muStart(mu), ng(ng_)
  { 
    
  }

  // Running in full SM nf=6
  AlphaQEDQCD(OSinput oi_, long double aMZ = pdg2014::aMZ, long double asMZ = pdg2014::asMZ): aStart(aMZ), asStart(asMZ), oi(oi_)
  {
    ng = 3;
    muStart = oi.MZ();
  }

  
  std::pair<long double,long double>  operator()(long double mu)
  {
    return runQEDQCD(aStart, asStart, muStart, mu, ng);    
  }

  long double QED(long double mu)
  {
    return runQEDQCD(aStart, asStart, muStart, mu, ng).second;
  }

  long double QCD(long double mu)
  {
    return runQEDQCD(aStart, asStart, muStart, mu, ng).first;
  }

};

#endif  // __BETAQEDQCD_HPP_
