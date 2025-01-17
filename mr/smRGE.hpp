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

#ifndef __BETASM_HPP__
#define __BETASM_HPP__

#include <iomanip>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>

#include "boost/numeric/odeint/integrate/integrate_adaptive.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/numeric/odeint/stepper/controlled_runge_kutta.hpp"

#include "tdecl.hpp"
#include "logger.hpp"
#include "constants.hpp"
#include "p2ms.hpp"

namespace mr
{
  using namespace boost::numeric::odeint;



  typedef runge_kutta_cash_karp54< SMCouplings,SMCouplings::value_type
                                   > error_stepper_type;
  typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;


  typedef std::vector<size_t> index_t;

  struct index_cmp_t : std::function<bool(index_t, index_t)> {
    bool operator ()(index_t const& a, index_t const& b) const {
      for (index_t::size_type i = 0; i < a.size(); ++i)
        if (a[i] != b[i])
          return a[i] < b[i];
      return false;
    }
  };




  class BetaSMFull
  {

    int NG;

    bool MultiplyByMinus1;
  
    std::map<index_t, Rt, index_cmp_t> be1;
    std::map<index_t, Rt, index_cmp_t> be2;
    std::map<index_t, Rt, index_cmp_t> be3;
    std::map<index_t, Rt, index_cmp_t> be4;
    std::map<index_t, Rt, index_cmp_t> be5;
    std::map<index_t, Rt, index_cmp_t> be6;
    std::map<index_t, Rt, index_cmp_t> be7;
  
    int pocoa1, pocoa2, pocoa3, pocoa4, pocoa5, pocoa6, pocoa7;
    size_t maxPower;

    void add(std::map<index_t, Rt, index_cmp_t>&, size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t, Rt);

  public:
  
    BetaSMFull(int pocoa1_, int pocoa2_, int pocoas_,
               int pocoat_, int pocoab_, int pocoatau_,
               int pocolam_,
               int NG_,
               bool MultiplyByMinus1_ = false);

    void operator() (const SMCouplings &, SMCouplings &, SMCouplings::value_type);

    // Forward direction, dir=False
    void setEvolutionDirection(bool dir)
    {
      MultiplyByMinus1 = dir;
    }

    Rt betaQCD(Rt);
  };









  // New class including mu2 and vev running 
  class BetaSM
  {

    std::unique_ptr<BetaSMFull> bSM;
    size_t ng;
  
    size_t maxPower;
    bool MultiplyByMinus1;
  
    int pocoa1, pocoa2, pocoa3, pocoa4, pocoa5, pocoa6, pocoa7, pocoa8, pocoa9;
    
    void add(std::map<index_t, Rt, index_cmp_t>&, size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t, Rt);

  public:
  
    BetaSM(int pocoa1_, int pocoa2_, int pocoas_,
           int pocoat_, int pocoab_, int pocoatau_,
           int pocolam_, int pocomu2_, int pocovev_,
           size_t NG_ = 3, bool MultiplyByMinus1_ = false);

    void operator() (const SMCouplings &, SMCouplings &, const Rt);
  
    void setEvolutionDirection(bool dir)
    {
      MultiplyByMinus1 = dir;
      bSM->setEvolutionDirection(dir);
    }

    static SMCouplings::value_type gamm2(const SMCouplings &a, size_t NG, int poco = 3);
    static SMCouplings::value_type gamv(const SMCouplings &a, size_t NG, int poco = 3);

  };


  // 
  // Couplings class
  // 
  template < int pocoa1, int pocoa2, int pocoas, int pocoat, int pocoab, int pocoatau, int pocolam > 
  class CouplingsSM
  {
    
    Rt                            mu0;
    SMCouplings                  aSM0;
    size_t                         NG;
    std::unique_ptr<BetaSMFull>   bep;
    
  public:
  
    CouplingsSM(Rt a1, Rt a2, Rt as, Rt at, Rt ab, Rt atau, Rt lam, Rt mu0_, size_t NG_) : mu0(mu0_), NG(NG_)
    {

      aSM0.push_back(a1);
      aSM0.push_back(a2);
      aSM0.push_back(as);
      aSM0.push_back(at);
      aSM0.push_back(ab);
      aSM0.push_back(atau);
      aSM0.push_back(lam);
      
      const size_t fw = 20;
      std::cout << std::scientific << std::setprecision(fw-10);
      std::cout << "# [Initial values:]" << std::endl;
      std::cout << "\t  a1 = " << std::setw(fw) << a1 
                << "   a2 = " << std::setw(fw) << a2 
                << "   as = " << std::setw(fw) << as <<std::endl; 
      std::cout << "\t  at = " << std::setw(fw) << at 
                << "   ab = " << std::setw(fw) << ab 
                << " atau = " << std::setw(fw) << atau <<std::endl; 
      std::cout << "\talam = " << std::setw(fw) << lam
                << "   MU = " << std::setw(fw) << sqrt(mu0)
                << "   NG = " << std::setw(fw) << NG <<std::endl; 
      
      bep = std::unique_ptr<BetaSMFull>(new  BetaSMFull(pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam, NG_));
    }


    SMCouplings operator()(SMCouplings::value_type mu2)
    {
        
      SMCouplings::value_type lEnd = log(mu2/mu0);

      bep->setEvolutionDirection(lEnd < 0);

      // Integration parameters
      // 
      // For the Runge-Kutta controller the error made during one step
      // is compared with 
      //    eps_abs + eps_rel * ( ax * |x| + adxdt * dt * |dxdt| ). 
      // If the error is smaller than this value the current
      // step is accepted, otherwise it is rejected and the step size is decreased.
      SMCouplings::value_type abs_err = 1.0e-22 , rel_err = 1.0e-20 , a_x = 1.0 , a_dxdt = 1.0;
      SMCouplings aSM(aSM0);    
    
      controlled_stepper_type 
        controlled_stepper(default_error_checker< SMCouplings::value_type, range_algebra, default_operations>
                           ( abs_err , rel_err , a_x , a_dxdt ) );
    
      integrate_adaptive( controlled_stepper , // Stepper function
                          *bep,                // Derivatives
                          aSM,                 // Initial values
                          Rt(0.0),             // t0 = Log[mu0/mu0]
                          Rt(std::abs(lEnd)),  // t  = Log[mu/mu0]
                          Rt(0.0001)           // Initial step size
                          );
      return aSM; 
    }
  
  
  };


  // 
  // Couplings and parameters
  // including [m] and [vev]
  template < int pocoa1, int pocoa2, int pocoas, int pocoat, int pocoab, int pocoatau, int pocolam, int pocomu2, int pocovev > 
  class ParametersSM
  {
  
    Rt                      mu0;
    SMCouplings            aSM0;
    size_t                   NG;
    std::unique_ptr<BetaSM> bep;
  
  public:
    ParametersSM(Rt a1, Rt a2, Rt as, Rt at, Rt ab, Rt atau, Rt lam, Rt mphi, Rt vev, Rt mu0_, size_t NG_ = 3) : mu0(mu0_), NG(NG_)
    {
    
      aSM0.push_back(a1);
      aSM0.push_back(a2);
      aSM0.push_back(as);
      aSM0.push_back(at);
      aSM0.push_back(ab);
      aSM0.push_back(atau);
      aSM0.push_back(lam);
      aSM0.push_back(mphi);
      aSM0.push_back(vev);
    
      const size_t fw = 20;
      lout(logINFO) << std::scientific << std::setprecision(fw-10);
      lout(logINFO) << "# [Initial values:]";
      lout(logINFO) << "\t  a1 = " << std::setw(fw) << a1 
                    << "   a2 = " << std::setw(fw) << a2 
                    << "   as = " << std::setw(fw) << as;
      lout(logINFO) << "\t  at = " << std::setw(fw) << at 
                    << "   ab = " << std::setw(fw) << ab 
                    << " atau = " << std::setw(fw) << atau;
      lout(logINFO) << "\talam = " << std::setw(fw) << lam
                    << "  vev = " << std::setw(fw) << vev
                    << " mphi = " << std::setw(fw) << mphi;
      lout(logINFO) << "\t  MU = " << std::setw(fw) << sqrt(mu0)
                    << "   NG = " << std::setw(fw) <<  NG; 
    
      bep = std::unique_ptr<BetaSM>(new BetaSM(pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam, pocomu2, pocovev, NG_));
    }

    // Constructor from Pole mass input
    template<class T>
    ParametersSM(const P2MS<T>& pi, size_t NG_ = 3) : NG(NG_)
    {

      mu0 = pow(pi.scale(),2);

      aSM0 = pi.ai();
    
      const size_t fw = 20;
      lout(logINFO) << std::scientific << std::setprecision(fw-10);
      lout(logINFO) << "# [Initial values:]";
      lout(logINFO) << "\t  a1 = " << std::setw(fw) << aSM0[0] 
                    << "   a2 = " << std::setw(fw) << aSM0[1]  
                    << "   as = " << std::setw(fw) << aSM0[2];
      lout(logINFO) << "\t  at = " << std::setw(fw) << aSM0[3]  
                    << "   ab = " << std::setw(fw) << aSM0[4]  
                    << " atau = " << std::setw(fw) << aSM0[5];
      lout(logINFO) << "\talam = " << std::setw(fw) << aSM0[6] 
                    << "  vev = " << std::setw(fw) << aSM0[8] 
                    << " mphi = " << std::setw(fw) << aSM0[7];
      lout(logINFO) << "\t  MU = " << std::setw(fw) << sqrt(mu0)
                    << "   NG = " << std::setw(fw) << NG;
    
      bep = std::unique_ptr<BetaSM>(new  BetaSM(pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam, pocomu2, pocovev, NG_));
    }


  
    SMCouplings operator()(const Rt& mu2End)
    {
    
      SMCouplings::value_type lEnd = log(mu2End/mu0);

      bep->setEvolutionDirection(lEnd < 0);
    
      // Integration parameters
      // 
      // For the Runge-Kutta controller the error made during one step
      // is compared with 
      //    eps_abs + eps_rel * ( ax * |x| + adxdt * dt * |dxdt| ). 
      // If the error is smaller than this value the current
      // step is accepted, otherwise it is rejected and the step size is decreased.
      SMCouplings::value_type  abs_err = 1.0e-12 , rel_err = 1.0e-5 , a_x = 1.0 , a_dxdt = 1.0;
      SMCouplings aSM(aSM0);

      lout(logDEBUG) << "                                             ";
      lout(logDEBUG) << "couplings values passed to the beta-functions";
      lout(logDEBUG) << "Initial scale mu0=" << sqrt(mu0);
      lout(logDEBUG) << "a1=" << aSM[0];
      lout(logDEBUG) << "a2=" << aSM[1];
      lout(logDEBUG) << "a3=" << aSM[2];
      lout(logDEBUG) << "at=" << aSM[3];
      lout(logDEBUG) << "ab=" << aSM[4];
      lout(logDEBUG) << "al=" << aSM[6];
      lout(logDEBUG) << "mu=" << aSM[7];
      lout(logDEBUG) << "                                             ";

      SMCouplings bSM(aSM.size());
      bep->operator()(aSM, bSM, 0);

      lout(logDEBUG) << "Starting values for Beta                     ";
      lout(logDEBUG) << "b1=" << bSM[0];
      lout(logDEBUG) << "b2=" << bSM[1];
      lout(logDEBUG) << "b3=" << bSM[2];
      lout(logDEBUG) << "bt=" << bSM[3];
      lout(logDEBUG) << "bb=" << bSM[4];
      lout(logDEBUG) << "bl=" << bSM[6];
    
      controlled_stepper_type 
        controlled_stepper(default_error_checker< SMCouplings::value_type , range_algebra , default_operations >
                           ( abs_err , rel_err , a_x , a_dxdt ) );
    
      integrate_adaptive( controlled_stepper , // Stepper function
                          std::ref(*bep),                // Derivatives
                          aSM,                 // Initial values
                          Rt(0.0),             // t0 = Log[mu0/mu0]
                          Rt(std::abs(lEnd)),  // t  = Log[mu/mu0]
                          Rt(0.001)            // Initial step size
                          );
    
      return aSM; 
    }

    // return couplings and its beta-functions
    std::pair<SMCouplings, SMCouplings> AandB(const Rt& mu2End)
    {
    
      Rt lEnd = log(mu2End/mu0);

      bep->setEvolutionDirection(lEnd < 0);

      // Integration parameters
      // 
      // For the Runge-Kutta controller the error made during one step
      // is compared with 
      //    eps_abs + eps_rel * ( ax * |x| + adxdt * dt * |dxdt| ). 
      // If the error is smaller than this value the current
      // step is accepted, otherwise it is rejected and the step size is decreased.
      SMCouplings::value_type abs_err = 1.0e-12 , rel_err = 1.0e-5 , a_x = 1.0 , a_dxdt = 1.0;
      SMCouplings aSM(aSM0);

      controlled_stepper_type 
        controlled_stepper(default_error_checker< SMCouplings::value_type , range_algebra , default_operations >
                           ( abs_err , rel_err , a_x , a_dxdt ) );
    
      integrate_adaptive( controlled_stepper , // Stepper function
                          *bep,                // Derivatives
                          aSM,                 // Initial values
                          Rt(0.0),             // t0 = Log[mu0/mu0]
                          Rt(fabs(lEnd)),      // t  = Log[mu/mu0]
                          Rt(0.001)            // Initial step size
                          );

      // beta-functions
      SMCouplings bSM(aSM0.size());

      bep->operator()(aSM,bSM,0);
    
      return std::make_pair(aSM,bSM); 
    }

  
  };
} // namespace mr


#endif  // __BETASM_HPP__

