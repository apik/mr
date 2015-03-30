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
// #include <boost/numeric/odeint/integrator_adaptive_stepsize.hpp>
// #include <boost/numeric/odeint.hpp>
#include "constants.hpp"
#include "p2ms.hpp"

using namespace boost::numeric::odeint;
typedef std::vector<long double> state_type;
typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

typedef std::vector<size_t> index_t;

struct index_cmp_t : std::binary_function<index_t, index_t, bool> {
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
  // Maximal numer of loops available for beta-functions
  static const size_t maxBetaOrder = 3;

  static const size_t maxPoco = maxBetaOrder + 2;

  bool MultiplyByMinus1;
  
  std::map<index_t, long double, index_cmp_t> be1;
  std::map<index_t, long double, index_cmp_t> be2;
  std::map<index_t, long double, index_cmp_t> be3;
  std::map<index_t, long double, index_cmp_t> be4;
  std::map<index_t, long double, index_cmp_t> be5;
  std::map<index_t, long double, index_cmp_t> be6;
  std::map<index_t, long double, index_cmp_t> be7;
  
  // int pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam;
  int pocoa1, pocoa2, pocoa3, pocoa4, pocoa5, pocoa6, pocoa7;
  size_t maxPower;

  void add(std::map<index_t, long double, index_cmp_t>&, size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t, long double);

public:
  
  BetaSMFull(int pocoa1_, int pocoa2_, int pocoas_,
             int pocoat_, int pocoab_, int pocoatau_,
             int pocolam_,
             int NG_,
             bool MultiplyByMinus1_ = false);

  void operator() (const state_type &, state_type &, const double);
  void multiplyByMinus1()
  {
    MultiplyByMinus1 = true;
  }
  long double betaQCD(long double);
};









// New class including mu2 and vev running 
class BetaSM
{

  BetaSMFull* bSM;
  size_t ng;
  
  size_t maxPower;
  bool MultiplyByMinus1;
  
  // int pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam;
  int pocoa1, pocoa2, pocoa3, pocoa4, pocoa5, pocoa6, pocoa7, pocoa8, pocoa9;

  void add(std::map<index_t, long double, index_cmp_t>&, size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t, long double);

public:
  
  BetaSM(int pocoa1_, int pocoa2_, int pocoas_,
         int pocoat_, int pocoab_, int pocoatau_,
         int pocolam_, int pocomu2_, int pocovev_,
         size_t NG_ = 3, bool MultiplyByMinus1_ = false);

  void operator() (const state_type &, state_type &, const double);
  
  void multiplyByMinus1()
  {
    MultiplyByMinus1 = true;
    bSM->multiplyByMinus1();
  }

  static long double bmu2(const state_type &a, size_t NG, int poco = 3);
  static long double gamv(const state_type &a, size_t NG, int poco = 3);

};



// couplings beta-functions and v.e.v.
class BetaVEV
{

  BetaSMFull* bSM;

  size_t ng;
  bool MultiplyByMinus1;
  
public:
  BetaVEV(size_t NG_, bool MultiplyByMinus1_ = false);

  void operator() (const state_type &, state_type &, const double);

  void multiplyByMinus1()
  {
    MultiplyByMinus1 = true;
    bSM->multiplyByMinus1();
  }

};


// couplings beta-functions and mu^2 
class BetaMu2
{

  BetaSMFull* bSM;

  size_t ng;
  bool MultiplyByMinus1;
  
public:
  BetaMu2(size_t NG_, bool MultiplyByMinus1_ = false);

  // static long double bmu2(const state_type &, size_t);
  void operator() (const state_type &, state_type &, const double);

  void multiplyByMinus1()
  {
    MultiplyByMinus1 = true;
    bSM->multiplyByMinus1();
  }

};


// couplings beta-functions, v.e.v. and mu^2 
class BetaVM
{

  BetaSMFull* bSM;
  size_t ng;
  
  size_t maxPower;
  bool MultiplyByMinus1;
  
public:
  BetaVM(size_t NG_, bool MultiplyByMinus1_ = false);

  void operator() (const state_type &, state_type &, const double);

  void multiplyByMinus1()
  {
    MultiplyByMinus1 = true;
    bSM->multiplyByMinus1();
  }

};







// 
// Couplings class
// 
template < int pocoa1, int pocoa2, int pocoas, int pocoat, int pocoab, int pocoatau, int pocolam > 
class CouplingsSM
{
  
  double      mu0;
  state_type aSM0;
  size_t       NG;
  BetaSMFull*   bep;

public:
  
  CouplingsSM(double a1, double a2, double as, double at, double ab, double atau, double lam, double mu0_, size_t NG_) : mu0(mu0_), NG(NG_)
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
      
      bep = new  BetaSMFull(pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam, NG_);
    }


  state_type operator()(long double mu2)
  {
    
    // BetaSMFull be = BetaSMFull(pocoa1, pocoa2, pocoas, pocoat, pocoab,
    //                 pocoatau, pocolam, NG);
    
    double lEnd = log(mu2/mu0);

    if (lEnd < 0) bep->multiplyByMinus1();
    // Integration parameters
    // 
    // For the Runge-Kutta controller the error made during one step
    // is compared with 
    //    eps_abs + eps_rel * ( ax * |x| + adxdt * dt * |dxdt| ). 
    // If the error is smaller than this value the current
    // step is accepted, otherwise it is rejected and the step size is decreased.
    double abs_err = 1.0e-12 , rel_err = 1.0e-10 , a_x = 1.0 , a_dxdt = 1.0;
    state_type aSM(aSM0);    

    controlled_stepper_type 
      controlled_stepper(default_error_checker< double , range_algebra , default_operations >
                         ( abs_err , rel_err , a_x , a_dxdt ) );
    
    integrate_adaptive( controlled_stepper , // Stepper function
                        *bep,            // Derivatives
                        aSM,             // Initial values
                        0.0,             // t0 = Log[mu0/mu0]
                        fabs(lEnd),      // t  = Log[mu/mu0]
                        0.0001           // Initial step size
                        );
    return aSM; 
  }
  
  
};


class CouplingsVevMu
{
  
  double      mu0;
  state_type aSM0;
  size_t       NG;
  BetaVM*     bep;
  
public:
  CouplingsVevMu(double a1, double a2, double as, double at, double ab, double atau, double lam, double vev, double mu2, double mu0_, size_t NG_ = 3) : mu0(mu0_), NG(NG_)
  {
    
    aSM0.push_back(a1);
    aSM0.push_back(a2);
    aSM0.push_back(as);
    aSM0.push_back(at);
    aSM0.push_back(ab);
    aSM0.push_back(atau);
    aSM0.push_back(lam);
    aSM0.push_back(vev);
    aSM0.push_back(mu2);
    
    const size_t fw = 20;
    std::cout << std::scientific << std::setprecision(fw-10);
    std::cout << "# [Initial values:]" << std::endl;
    std::cout << "\t  a1 = " << std::setw(fw) << a1 
              << "   a2 = " << std::setw(fw) << a1 
              << "   as = " << std::setw(fw) << as <<std::endl; 
    std::cout << "\t  at = " << std::setw(fw) << at 
              << "   ab = " << std::setw(fw) << ab 
              << " atau = " << std::setw(fw) << atau <<std::endl; 
    std::cout << "\talam = " << std::setw(fw) << lam
              << "  vev = " << std::setw(fw) << vev
              << "  mu2 = " << std::setw(fw) << mu2 <<std::endl; 
    std::cout << "\t  MU = " << std::setw(fw) << sqrt(mu0)
              << "   NG = " << std::setw(fw) << NG <<std::endl; 
    
    bep = new  BetaVM(NG_);
  }
  
  
  state_type operator()(const long double& mu2End)
  {
    
    // BetaSMFull be = BetaSMFull(pocoa1, pocoa2, pocoas, pocoat, pocoab,
    //                 pocoatau, pocolam, NG);
    
    double lEnd = log(mu2End/mu0);

    if (lEnd < 0) bep->multiplyByMinus1();
    
    // Integration parameters
    // 
    // For the Runge-Kutta controller the error made during one step
    // is compared with 
    //    eps_abs + eps_rel * ( ax * |x| + adxdt * dt * |dxdt| ). 
    // If the error is smaller than this value the current
    // step is accepted, otherwise it is rejected and the step size is decreased.
    double abs_err = 1.0e-12 , rel_err = 1.0e-10 , a_x = 1.0 , a_dxdt = 1.0;
    state_type aSM(aSM0);    
    // controlled_stepper_standard< stepper_rk5_ck< state_type > >
    //   controlled_rk5( abs_err , rel_err , a_x , a_dxdt );
 

    // integrate_adaptive( controlled_rk5 , // Stepper function
    //                     *bep,            // Derivatives
    //                     aSM,             // Initial values
    //                     0.0,             // t0 = Log[mu0/mu0]
    //                     lEnd,            // t  = Log[mu/mu0]
    //                     0.0001             // Initial step size
    //                     );

    controlled_stepper_type 
      controlled_stepper(default_error_checker< double , range_algebra , default_operations >
                         ( abs_err , rel_err , a_x , a_dxdt ) );
    
    integrate_adaptive( controlled_stepper , // Stepper function
                        *bep,            // Derivatives
                        aSM,             // Initial values
                        0.0,             // t0 = Log[mu0/mu0]
                        fabs(lEnd),            // t  = Log[mu/mu0]
                        0.0001             // Initial step size
                        );
    
    return aSM; 
  }
  
  
};

// template < int pocoa1, int pocoa2, int pocoas, int pocoat, int pocoab, int pocoatau, int pocolam, int pocomu2> 
class CouplingsMu
{
  
  double      mu0;
  state_type aSM0;
  size_t       NG;
  BetaMu2*     bep;

public:
  CouplingsMu(double a1, double a2, double as, double at, double ab, double atau, double lam, double mu2, double mu0_, size_t NG_ = 3) : mu0(mu0_), NG(NG_)
    {

      aSM0.push_back(a1);
      aSM0.push_back(a2);
      aSM0.push_back(as);
      aSM0.push_back(at);
      aSM0.push_back(ab);
      aSM0.push_back(atau);
      aSM0.push_back(lam);

      aSM0.push_back(mu2);
      
      const size_t fw = 20;
      std::cout << std::scientific << std::setprecision(fw-10);
      std::cout << "# [Initial values:]" << std::endl;
      std::cout << "\t  a1 = " << std::setw(fw) << a1 
                << "   a2 = " << std::setw(fw) << a1 
                << "   as = " << std::setw(fw) << as <<std::endl; 
      std::cout << "\t  at = " << std::setw(fw) << at 
                << "   ab = " << std::setw(fw) << ab 
                << " atau = " << std::setw(fw) << atau <<std::endl; 
      std::cout << "\talam = " << std::setw(fw) << lam
      
                << "  mu2 = " << std::setw(fw) << mu2 <<std::endl; 
      std::cout << "\t  MU = " << std::setw(fw) << sqrt(mu0)
                << "   NG = " << std::setw(fw) << NG <<std::endl; 

       bep = new  BetaMu2(NG_);
    }

  state_type operator()(const long double& mu2End)
  {
    
    // BetaSMFull be = BetaSMFull(pocoa1, pocoa2, pocoas, pocoat, pocoab,
    //                 pocoatau, pocolam, NG);

    double lEnd = log(mu2End/mu0);

    if (lEnd < 0) bep->multiplyByMinus1();
    
    // Integration parameters
    // 
    // For the Runge-Kutta controller the error made during one step
    // is compared with 
    //    eps_abs + eps_rel * ( ax * |x| + adxdt * dt * |dxdt| ). 
    // If the error is smaller than this value the current
    // step is accepted, otherwise it is rejected and the step size is decreased.
    double abs_err = 1.0e-12 , rel_err = 1.0e-10 , a_x = 1.0 , a_dxdt = 1.0;
    state_type aSM(aSM0);    
    // controlled_stepper_standard< stepper_rk5_ck< state_type > >
    //   controlled_rk5( abs_err , rel_err , a_x , a_dxdt );
 

    // integrate_adaptive( controlled_rk5 , // Stepper function
    //                     *bep,            // Derivatives
    //                     aSM,             // Initial values
    //                     0.0,             // t0 = Log[mu0/mu0]
    //                     lEnd,            // t  = Log[mu/mu0]
    //                     0.0001             // Initial step size
    //                     );
    controlled_stepper_type 
      controlled_stepper(default_error_checker< double , range_algebra , default_operations >
                         ( abs_err , rel_err , a_x , a_dxdt ) );

    integrate_adaptive( controlled_stepper , // Stepper function
                        *bep,            // Derivatives
                        aSM,             // Initial values
                        0.0,             // t0 = Log[mu0/mu0]
                        fabs(lEnd),            // t  = Log[mu/mu0]
                        0.0001             // Initial step size
                        );

    return aSM; 
  }


};



template < int pocoa1, int pocoa2, int pocoas, int pocoat, int pocoab, int pocoatau, int pocolam, int pocomu2, int pocovev > 
class Couplings
{
  
  double      mu0;
  state_type aSM0;
  size_t       NG;
  BetaSM*     bep;
  
public:
  Couplings(double a1, double a2, double as, double at, double ab, double atau, double lam, double mu2, double vev, double mu0_, size_t NG_ = 3) : mu0(mu0_), NG(NG_)
  {
    
    aSM0.push_back(a1);
    aSM0.push_back(a2);
    aSM0.push_back(as);
    aSM0.push_back(at);
    aSM0.push_back(ab);
    aSM0.push_back(atau);
    aSM0.push_back(lam);
    aSM0.push_back(mu2);
    aSM0.push_back(vev);
    
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
              << "  vev = " << std::setw(fw) << vev
              << "  mu2 = " << std::setw(fw) << mu2 <<std::endl; 
    std::cout << "\t  MU = " << std::setw(fw) << sqrt(mu0)
              << "   NG = " << std::setw(fw) << NG <<std::endl; 
    
    bep = new  BetaSM(pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam, pocomu2, pocovev, NG_);
  }

  // Constructor from Pole mass input
  
  Couplings(const P2MS& pi, size_t NG_ = 3) : NG(NG_)
  {

    mu0 = pow(pi.scale(),2);

    aSM0 = pi.ai();
    // aSM0.push_back(a1);
    // aSM0.push_back(a2);
    // aSM0.push_back(as);
    // aSM0.push_back(at);
    // aSM0.push_back(ab);
    // aSM0.push_back(atau);
    // aSM0.push_back(lam);
    // aSM0.push_back(mu2);
    // aSM0.push_back(vev);
    
    const size_t fw = 20;
    std::cout << std::scientific << std::setprecision(fw-10);
    std::cout << "# [Initial values:]" << std::endl;
    std::cout << "\t  a1 = " << std::setw(fw) << aSM0[0] 
              << "   a2 = " << std::setw(fw) << aSM0[1]  
              << "   as = " << std::setw(fw) << aSM0[2]  <<std::endl; 
    std::cout << "\t  at = " << std::setw(fw) << aSM0[3]  
              << "   ab = " << std::setw(fw) << aSM0[4]  
              << " atau = " << std::setw(fw) << aSM0[5]  <<std::endl; 
    std::cout << "\talam = " << std::setw(fw) << aSM0[6] 
              << "  vev = " << std::setw(fw) << aSM0[8] 
              << "  mu2 = " << std::setw(fw) << aSM0[7]  <<std::endl; 
    std::cout << "\t  MU = " << std::setw(fw) << sqrt(mu0)
              << "   NG = " << std::setw(fw) << NG <<std::endl; 
    
    bep = new  BetaSM(pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam, pocomu2, pocovev, NG_);
  }


  
  state_type operator()(const long double& mu2End)
  {
    
    double lEnd = log(mu2End/mu0);

    if (lEnd < 0) bep->multiplyByMinus1();
    
    // Integration parameters
    // 
    // For the Runge-Kutta controller the error made during one step
    // is compared with 
    //    eps_abs + eps_rel * ( ax * |x| + adxdt * dt * |dxdt| ). 
    // If the error is smaller than this value the current
    // step is accepted, otherwise it is rejected and the step size is decreased.
    double abs_err = 1.0e-12 , rel_err = 1.0e-10 , a_x = 1.0 , a_dxdt = 1.0;
    state_type aSM(aSM0);    

    controlled_stepper_type 
      controlled_stepper(default_error_checker< double , range_algebra , default_operations >
                         ( abs_err , rel_err , a_x , a_dxdt ) );
    
    integrate_adaptive( controlled_stepper , // Stepper function
                        *bep,            // Derivatives
                        aSM,             // Initial values
                        0.0,             // t0 = Log[mu0/mu0]
                        fabs(lEnd),      // t  = Log[mu/mu0]
                        0.0001           // Initial step size
                        );
    
    return aSM; 
  }
  
  
};



#endif  // __BETASM_HPP__

