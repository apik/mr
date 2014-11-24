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

#ifndef __BETAQCD_HPP_
#define __BETAQCD_HPP_

#include <boost/numeric/odeint/integrator_adaptive_stepsize.hpp>

#include <stdexcept>

#include "sminput.hpp"
#include "constants.hpp"


using namespace boost::numeric::odeint;

typedef std::vector<long double> state_type;

class BetaQCD
{
  
  size_t loops;
  double nf;
  
  bool MultiplyByMinus1;
public:
  BetaQCD( size_t loops_, double nf_, bool MultiplyByMinus1_ = false) : loops(loops_), nf(nf_), MultiplyByMinus1(MultiplyByMinus1_){ }
  double b0()
  {
    return loops > 0 ? 
      (11.-2./3.*nf) : 0;
  }
  double b1()
  {
    return loops > 1 ? 
      (102.-38./3.*nf) : 0;
  }
  double b2()
  {
    return loops > 2 ? 
      (2857./2.-5033./18.*nf + 325./54.*nf*nf) : 0;
  }
  double b3()
  {
    return loops > 3 ? 
      (149753./6. + 3564.*Zeta3
      - (1078361./162. + 6508./27.*Zeta3)*nf
      + (50065./162. + 6472./81.*Zeta3)*nf*nf
       + 1093./729.*nf*nf*nf) : 0;
  }
  void operator() (const state_type &a, state_type &dadt, const double t)
  {
    double minusC = MultiplyByMinus1 ? -1. : 1.;
    dadt[0] =  minusC         // We need -\\beta if evolve from t=0 to t < 0
      *(-b0()*pow(a[0],2) - b1()*pow(a[0],3) - b2()*pow(a[0],4) - b3()*pow(a[0],5));
  }
};


// run asStart from muStart to muEnd, 
// using RGE with nf active flavours 
double run(long double asStart, long double muStart, long double muEnd, size_t NF = 5, size_t loops = 4);
//double run(long double asStart, long double muStart, long double muEnd, size_t NF, size_t loops);




/*
  
         ^
         |
  nf = 6 +                      o
         |                     /
  nf = 5 +     o---<---#--->---
         |
         |
         +-----+-------+--------+--->
        0      Mb      MZ       Mt
        
        We use as(MZ) as initial input, than 
        for applications in MR we need:
        1. as(Mb,nf=5) we use nf=5 running from MZ
        2. as(Mt,nf=6) and as(mu>Mt,nf=6) we use 
        nf=5 running from MZ to Mt and at threshold 
        mu=Mt match as(nf=5) and as(nf=6) than run 
        to higher scale with as(nf=6) if needed.
 */

long double as5nf2as6nf(long double M, long double mu, long double as, size_t nl = 5, size_t ord = 3);
// long double as5nf2as6nf(long double, long double, long double, size_t, size_t);



class AlphaS
{
  double muStart;
  double asStart;
  size_t   loops;
  size_t nfFixed;
  OSinput     oi;

public:

  // Running with fixed nf
  AlphaS(long double asMZ = pdg2014::asMZ, long double mu = pdg2014::MZ, size_t loops_ = 4, size_t nfFixed_ = 5) : asStart(asMZ), muStart(mu),loops(loops_), nfFixed(nfFixed_)
  { 
    
  }

  // Running down to bottom mass with nf=5
  // and upto Mt with threshold at Mt
  AlphaS(OSinput oi_, long double asMZ = pdg2014::asMZ, size_t loops_ = 4, size_t nfFixed_ = 0) 
    : loops(loops_), asStart(asMZ), oi(oi_), nfFixed(nfFixed_)
  { 
    muStart = oi.MZ();

    // nf is not fixed
    // nfFixed = 0, default
  }

  
  double operator()(long double mu)
  {
    
    // Running with decoupling
    if(nfFixed == 0)
      {

        // Run only down to Mb
        if( mu < oi.Mb() ) 
          throw std::logic_error("ERROR: running with nf = 5 to scale mu < Mb"); 
        
        // Run only up to Mt with nf=5
        else if( mu < oi.Mt() ) 
          {
            return run(asStart, muStart, mu, 5 );
          }
        
        // Threshold at Mt
        else 
          {
            long double asMt5 = run(asStart, muStart, oi.Mt(), 5 );
            
            // We use 3-loop decoupling for 4-loop running
            long double asMt6 = as5nf2as6nf(oi.Mt(), oi.Mt(), asMt5, /* nl= */5, 3);
            
            // Return as(nf=6,mu=Mt)
            if (mu == oi.Mt())
              return asMt6;

            // Run with nf=6 up to mu > Mt
            else
              return run(asMt6, oi.Mt(), mu, 6 );
              
          }
      }
    // Running with fixed nf
    else
      return run(asStart, muStart, mu, nfFixed);    
  }


  long double upto(long double mu2)
  {
    
  }

};

#endif  // __BETAQCD_HPP_
