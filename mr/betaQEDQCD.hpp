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

#include <boost/numeric/odeint/integrator_adaptive_stepsize.hpp>

#include <stdexcept>

#include "sminput.hpp"
#include "constants.hpp"
#include "betaQCD.hpp"

using namespace boost::numeric::odeint;

typedef std::vector<long double> state_type;

class BetaQEDQCD
{
  
  BetaQCD* pbQCD;
  
  size_t nu,nd,nl;
  
  bool MultiplyByMinus1;
public:
  BetaQEDQCD( size_t nu_,size_t nd_,size_t nl_, bool MultiplyByMinus1_ = false) : nu(nu_),nd(nd_),nl(nl_), MultiplyByMinus1(MultiplyByMinus1_){
    pbQCD = new BetaQCD(4, nu + nd, MultiplyByMinus1);
  }
  
  // double b3()
  // {
  //   return loops > 3 ? 
  //     (149753./6. + 3564.*Zeta3
  //     - (1078361./162. + 6508./27.*Zeta3)*nf
  //     + (50065./162. + 6472./81.*Zeta3)*nf*nf
  //      + 1093./729.*nf*nf*nf)/256. : 0;
  // }
  void operator() (const state_type &av4pi, state_type &dadt, const double t)
  {

    long double Nu = static_cast<long double>(nu);
    long double Nd = static_cast<long double>(nd);
    long double Nl = static_cast<long double>(nl);
    
    double minusC = MultiplyByMinus1 ? -1. : 1.;

    if(av4pi.size() != 2)
      throw std::logic_error("ERROR: for mu^2 running 2 constants in input needed: aQCD and aQED");

    pbQCD->operator()(av4pi, dadt, t);

    // std::cout << "Beta QCD:" << dadt[0]<< std::endl; 
    long double as4pi = av4pi[0], a4pi = av4pi[1];
    // QCD beta-function
    // dadt[0] =                   // We need -\\beta if evolve from t=0 to t < 0
    //   minusC*((pow(as4pi,2)) * (-11 + (2/3.) * (Nd) + (2/3.) * (Nu)) + 
    //           (pow(a4pi,2)) * ((pow(as4pi,2)) * ((-1/81.) * (Nd) + (-22/243.) * (pow(Nd,2)) + (-22/81.) * ((Nd) * (Nl))
    //                                              + (-16/81.) * (Nu) + (-176/243.) * ((Nd) * (Nu))
    //                                              + (-88/81.) * ((Nl) * (Nu)) + (-352/243.) * (pow(Nu,2)))) +
    //           (pow(as4pi,2)) * ((a4pi) * ((2/9.) * (Nd) + (8/9.) * (Nu)) +
    //                             (as4pi) * (-102 + (38/3.) * (Nd) + (38/3.) * (Nu))) + 
    //           (pow(as4pi,3)) * ((a4pi) * ((28/27.) * (Nd) + (112/27.) * (Nu)) + 
    //                             (as4pi) * (-2857/2. + (5033/18.) * (Nd) + (-325/54.) * (pow(Nd,2)) + (5033/18.) * (Nu)
    //                                        + (-325/27.) * ((Nd) * (Nu)) + (-325/54.) * (pow(Nu,2)))));


    
    // QED beta-function
    dadt[1] =                   // We need -\\beta if evolve from t=0 to t < 0
      minusC*((pow(a4pi,2)) * ((4/9.) * (Nd) + (4/3.) * (Nl) + (16/9.) * (Nu))
              // + 
              // (pow(a4pi,2)) * ((pow(as4pi,2)) * ((500/27.) * (Nd) + (-88/81.) * (pow(Nd,2))+ (2000/27.) * (Nu)
              //                                    + (-440/81.) * ((Nd) * (Nu)) +  (-352/81.) * (pow(Nu,2))))+
              // (pow(a4pi,2)) * ((a4pi) * ((4/27.) *  (Nd) + (4) * (Nl) + (64/27.) * (Nu)) +
              //                  (as4pi) * ((16/9.) * (Nd) + (64/9.) * (Nu))) +
              // (pow(a4pi,3)) * ((as4pi) * ((-16/81.) * (Nd) + (-256/81.) * (Nu)) +
              //                  (a4pi) * ((-2/243.) * (Nd) + (-44/729.) * (pow(Nd,2)) + (-2) * (Nl) + (-440/243.) * ((Nd) * (Nl))
              //                            + (-44/9.) * (pow(Nl,2)) + (-128/243.) * (Nu) + (-880/729.) * ((Nd) * (Nu))
              //                            +  (-2288/243.) * ((Nl) * (Nu)) + (-2816/729.) * (pow(Nu,2))))
              );
    
  }
};



std::pair<double,double> runQEDQCD(long double aStart, long double asStart, long double muStart, long double muEnd, size_t nu = 2, size_t nd = 3, size_t nl = 3);



class AlphaQEDQCD
{
  double muStart;
  double  aStart;
  double asStart;
  size_t   loops;
  size_t nuFixed;
  size_t ndFixed;
  size_t nlFixed;
  OSinput     oi;

public:

  // Running with fixed nf
  AlphaQEDQCD(long double aMZ = pdg2014::aMZ,long double asMZ = pdg2014::asMZ, long double mu = pdg2014::MZ, size_t nuFixed_ = 2, size_t ndFixed_ = 3, size_t nlFixed_ = 6) : asStart(asMZ), muStart(mu), nuFixed(nuFixed_), ndFixed(ndFixed_), nlFixed(nlFixed_)
  { 
    
  }

  // Running down to bottom mass with nf=5
  // and upto Mt with threshold at Mt
  AlphaQEDQCD(OSinput oi_, long double aMZ = pdg2014::aMZ, long double asMZ = pdg2014::asMZ, size_t nuFixed_ = 0, size_t ndFixed_ = 0, size_t nlFixed_ = 0) 
    : aStart(aMZ), asStart(asMZ), oi(oi_), nuFixed(nuFixed_), ndFixed(ndFixed_), nlFixed(nlFixed_)
  { 
    muStart = oi.MZ();

    // nf is not fixed
    // nfFixed = 0, default
  }

  
  double operator()(long double mu)
  {
    
    // Running with decoupling
    if(nuFixed == 0 &&
       ndFixed == 0 &&
       nlFixed == 0)
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
    // else
    //   return run(asStart, muStart, mu, nfFixed);    
  }

};

#endif  // __BETAQEDQCD_HPP_
