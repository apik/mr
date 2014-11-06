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
#include <boost/numeric/odeint.hpp>
#include "constants.hpp"

using namespace boost::numeric::odeint;

typedef std::vector<long double> state_type;

class BetaQCD
{

  size_t loops;
  double nf;

public:
  BetaQCD( size_t loops_, double nf_) : loops(loops_), nf(nf_) { }
  
  double b0()
  {
    return loops > 0 ? 
      (11.-2./3.*nf)/4. : 0;
  }
  double b1()
  {
    return loops > 1 ? 
      (102.-38./3.*nf)/16. : 0;
  }
  double b2()
  {
    return loops > 2 ? 
      (2857./2.-5033./18.*nf + 325./54.*nf*nf)/64. : 0;
  }
  double b3()
  {
    return loops > 3 ? 
      (149753./6. + 3564.*Zeta3
      - (1078361./162. + 6508./27.*Zeta3)*nf
      + (50065./162. + 6472./81.*Zeta3)*nf*nf
       + 1093./729.*nf*nf*nf)/256. : 0;
  }
  void operator() (const state_type &a, state_type &dadt, const double t)
  {
    dadt[0] = -b0()*pow(a[0],2) - b1()*pow(a[0],3) - b2()*pow(a[0],4) - b3()*pow(a[0],5);
  }
};


class AlphaS
{
  double   mu0;
  size_t loops;
  size_t    nf;
public:
  AlphaS(size_t loops_ = 4, size_t nf_ = 5) : loops(loops_), nf(nf_)
  { }
  double operator()(long double mu2)
  {
    state_type a0(1);

    // Starting value
    a0[0] = 0.1184/Pi; 
            
    BetaQCD beta4l5nf(loops,nf);
    
    double mu0 = pow(91.1876,2);
    
    double lEnd = log(mu2/mu0);
    
    controlled_stepper_standard< stepper_rk5_ck< state_type > >
      controlled_rk5( 1E-6 , 1E-7 , 1.0 , 1.0 );
 
    integrate_adaptive( controlled_rk5 ,
                        beta4l5nf, a0, 0.0, lEnd, 0.01  );
    
    return a0[0]*Pi; 

  }
};

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

long double as5nf2as6nf(long double MM, long double mu2, long double as4Pi, size_t nl = 5, size_t ord = 4)
{

  long double Lmu2MM = log(mu2/MM);

  // z[as(nf=5)] = as(nf=6)/as(nf=5)
  long double z = 1;

  if(ord > 0)

    z += as4Pi*Lmu2MM * ( 2./3. );

  if(ord > 1)
    {
      z +=  + pow(as4Pi,2) * ( 14./3. );
      
      z +=  + pow(as4Pi,2)*Lmu2MM * ( 38./3. );
      
      z +=  + pow(as4Pi,2)*pow(Lmu2MM,2) * ( 4./9. );
      
    }
  
  if(ord > 2)
    {
      
      z +=  + pow(as4Pi,3) * ( 58933./1944. + 80507./432.*Zeta3 + 128./3.*Zeta2 + 128./9.*log(2)*Zeta2 );
      
      z +=  + pow(as4Pi,3)*nl * (  - 2479./486. - 64./9.*Zeta2 );
      
      z +=  + pow(as4Pi,3)*Lmu2MM * ( 8941./27. );
      
      z +=  + pow(as4Pi,3)*Lmu2MM*nl * (  - 409./27. );
      
      z +=  + pow(as4Pi,3)*pow(Lmu2MM,2) * ( 511./9. );

      z +=  + pow(as4Pi,3)*pow(Lmu2MM,3) * ( 8./27. );
    }

  if(ord > 3)
    {

      z +=  + pow(as4Pi,4) * ( 592371712./382725. - 21814592./2835.*
                               a5 - 25433192./1701.*a4 + 40596749./5670.*Zeta5 + 71102219./8505.
                               *Zeta4 + 2408412383./340200.*Zeta3 + 11153936./1215.*Zeta2 - 46048./
                               27.*Zeta2*Zeta3 - 18636934./2835.*log(2)*Zeta4 - 131456./81.*log(2)*Zeta2 + 
                               5826074./1701.*pow(log(2),2)*Zeta2 - 5453648./8505.*pow(log(2),3)*Zeta2
                               - 3179149./5103.*pow(log(2),4) + 2726824./42525.*pow(log(2),5) );
      
      z +=  + pow(as4Pi,4)*nl * (  - 1773073./2916. - 692./81.*a4 - 
                                   460./9.*Zeta5 + 697709./648.*Zeta4 - 4756441./3888.*Zeta3 - 71296./
                                   81.*Zeta2 - 5632./81.*log(2)*Zeta2 + 1709./81.*pow(log(2),2)*Zeta2 - 173./
                                   486.*pow(log(2),4) );
      
      z +=  + pow(as4Pi,4)*pow(nl,2) * ( 140825./5832. + 76./27.*Zeta3
                                         + 1664./81.*Zeta2 );
      
      z +=  + pow(as4Pi,4)*Lmu2MM * ( 21084715./2916. + 2922161./648.*Zeta3
                                 + 8960./9.*Zeta2 + 8960./27.*log(2)*Zeta2 );
      
      z +=  + pow(as4Pi,4)*Lmu2MM*nl * (  - 1140191./1458. - 132283./324.*
                                     Zeta3 - 6016./27.*Zeta2 - 512./27.*log(2)*Zeta2 );
      
      z +=  + pow(as4Pi,4)*Lmu2MM*pow(nl,2) * ( 1679./729. + 256./27.*Zeta2);
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,2) * ( 94078./27. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,2)*nl * (  - 18230./81. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,2)*pow(nl,2) * ( 493./81. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,3) * ( 28298./81. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,3)*nl * (  - 428./27. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,4) * ( 16./81. );
      
    }

  if(ord > 4)
    std::cerr << "WARNING: Only 4-loop decoupling relations available."<< std::endl;

  return as4Pi*z;
}

#endif  // __BETAQCD_HPP_
