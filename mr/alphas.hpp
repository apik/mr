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

                       
#include <stdexcept>
#include <vector>
#include "sminput.hpp"
#include "constants.hpp"
#include "tdecl.hpp"


namespace mr
{
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
    double b4()
    {
      return loops > 4 ? 
        (1205./2916. - 152./81.*Zeta3)*pow(nf,4) +
        (-630559./5832. - 48722./243.*Zeta3 + 1618./27.*Zeta4 + 460./9.*Zeta5)*pow(nf,3) + 
        (25960913./1944. + 698531./81.*Zeta3 - 10526./9.*Zeta4 - 381760./81.*Zeta5)*pow(nf,2) +
        (-336460813./1944. - 4811164./81.*Zeta3 + 33935./6.*Zeta4 + 1358995./27.*Zeta5)*pow(nf,1) +
        (8157455./16. + 621885./2.*Zeta3 - 88209./2.*Zeta4 - 288090.*Zeta5)*pow(nf,0) : 0;
    }

    void operator() (const SMCouplings &a, SMCouplings &dadt, const double t)
    {
      double minusC = MultiplyByMinus1 ? -1. : 1.;
      dadt[0] =  minusC         // We need -\\beta if evolve from t=0 to t < 0
        *(-b0()*pow(a[0],2) - b1()*pow(a[0],3) - b2()*pow(a[0],4) - b3()*pow(a[0],5) - b4()*pow(a[0],6));
    }
  };


  // run asStart from muStart to muEnd, 
  // using RGE with nf active flavours 
  double run(double asStart, double muStart, double muEnd, size_t NF = 5, size_t loops = 4);
  //double run(double asStart, double muStart, double muEnd, size_t NF, size_t loops);




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

  double as5nf2as6nf(double M, double mu, double as, size_t nl = 5, size_t ord = 3);
  // double as5nf2as6nf(double, double, double, size_t, size_t);



  class AlphaS
  {
    double muStart;
    double asStart;
    size_t   loops;
    size_t nfFixed;
    OSinput     oi;
    // Thresholds
    double    mbth;
    double    mtth;

  public:

    // Running down to bottom mass with nf=5
    // and upto Mt with threshold at Mt
    AlphaS(OSinput oi_, double asMZ = pdg2014::asMZ, size_t loops_ = 4, size_t nfFixed_ = 0) 
      : loops(loops_), asStart(asMZ), oi(oi_), nfFixed(nfFixed_)
    { 
      muStart = oi.MZ();
      mbth    = oi.Mb();
      mtth    = oi.Mt();
      // nf is not fixed
      // nfFixed = 0, default
    }

    // Same with manual thresholds for Mb and Mt
    AlphaS(double mu = pdg2014::MZ, double asMZ = pdg2014::asMZ, size_t loops_ = 4, size_t nfFixed_ = 5, double mtth_ = pdg2014::Mt) 
      : muStart(mu), mtth(mtth_), loops(loops_), asStart(asMZ),  nfFixed(nfFixed_)
    { 
      mbth = 0;
    }

  
    double operator()(double mu)
    {
    
      // Running with decoupling
      if(nfFixed == 0)
        {

          // Run only down to Mb
          if( mu < mbth ) 
            throw std::logic_error("ERROR: running with nf = 5 to scale mu < Mb"); 
        
          // Run only up to Mt with nf=5
          else if( mu < mtth ) 
            {
              return run(asStart, muStart, mu, 5, loops );
            }
        
          // Threshold at Mt
          else 
            {
              double asMt5 = run(asStart, muStart, mtth, 5, loops );
            
              // We use (L-1)-loop decoupling for L-loop running
              double asMt6 =
                as5nf2as6nf(mtth, mtth, asMt5, /* nl= */5, loops - 1);
            
              // Return as(nf=6,mu=Mt)
              if (mu == mtth)
                return asMt6;

              // Run with nf=6 up to mu > Mt
              else
                return run(asMt6, mtth, mu, 6, loops );
              
            }
        }
      // Running with fixed nf
      else
        return run(asStart, muStart, mu, nfFixed, loops);    
    }

  };

} // namespace mr

#endif  // __BETAQCD_HPP_
