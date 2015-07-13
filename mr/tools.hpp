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

#ifndef __TOOLS_HPP__
#define __TOOLS_HPP__

#include "sminput.hpp"
#include "betaQCD.hpp"
#include "betaSM.hpp"

#include <Eigen/Dense>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }

};


// Meta-stability - Stability bound
struct Mt_Stability : Functor<double>
{
  OSinput oi;
  long double aSMZ;
  
  Mt_Stability(void): Functor<double>(2,2) {}

  // Explicit dependence on asMZ is usefull for plots 
  // depending on asMZ errors
  Mt_Stability(const OSinput& in_, long double alphaS_ = pdg2014::asMZ) : oi(in_), aSMZ(alphaS_), Functor<double>(2,2)
  {
  }

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {
    OSinput oiMt(oi);
    oiMt.setMt(x(0));

    AlphaS as(oiMt,aSMZ);

    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMt,pdg2014::Gf, as(oiMt.Mt()), oiMt.Mt(), order::all);
    
    Couplings<3,3,3,
      3,3,-1,
      3,0,0> avP2MS(pMSmt);
    long double mu2 = pow(10.,x(1));
    
    lout(logDEBUG) << "Instability scale Lambda_I = " << sqrt(mu2);


    std::pair<state_type, state_type> ab = avP2MS.AandB(mu2);
    fvec(0) = pow(ab.first[couplings::lam],1);
    fvec(1) = pow(ab.second[couplings::lam],1);

    lout(logDEBUG) << "lam  = " << ab.first[couplings::lam];
    lout(logDEBUG) << "blam = " << ab.second[couplings::lam];
    return 0;
  }
};




// Criterium for Meta-stability-Instability bound
struct Mt_Instability : Functor<double>
{
  OSinput oi;
  long double aSMZ;
  
  Mt_Instability(void): Functor<double>(2,2) {}
  
  // Explicit dependence on asMZ is usefull for plots 
  // depending on asMZ errors
  Mt_Instability(const OSinput& in_, long double alphaS_ = pdg2014::asMZ) : oi(in_), aSMZ(alphaS_), Functor<double>(2,2)
  {
  }

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {
    OSinput oiMt(oi);
    oiMt.setMt(x(0));

    AlphaS as(oiMt, aSMZ);

    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMt,pdg2014::Gf,  as(oiMt.Mt()), oiMt.Mt(), order::all);
    
    Couplings <3,3,3,
      3,3,-1,
      3,0,0> avP2MS(pMSmt);
    long double mu2 = pow(10.,x(1));
    
    lout(logDEBUG) << "Instability scale Lambda_I = " << sqrt(mu2);


    std::pair<state_type, state_type> ab = avP2MS.AandB(mu2);

    // eq.from 1408.0292
    long double LamMin = 1./(-14.53 - 0.153*log(sqrt(mu2)));
    lout(logINFO) << "LamMin  = " << LamMin;

    fvec(0) = ab.first[couplings::lam] - LamMin;  // lambda=lambda_min
    fvec(1) = ab.second[couplings::lam]; // beta = 0

    lout(logDEBUG) << "lam  = " << ab.first[couplings::lam];
    lout(logDEBUG) << "blam = " << ab.second[couplings::lam];
    return 0;
  }
};


// Meta-stability - Stability bound
struct Mt_Lambda0 : Functor<double>
{
  OSinput oi;
  long double mu;
  long double aSMZ;
  Mt_Lambda0(void): Functor<double>(1,1) {}

  // Explicit dependence on asMZ is usefull for plots 
  // depending on asMZ errors
  Mt_Lambda0(const OSinput& in_, long double mu_, long double alphaS_ = pdg2014::asMZ) : oi(in_), mu(mu_), aSMZ(alphaS_), Functor<double>(1,1)
  {
  }

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {
    OSinput oiMt(oi);
    oiMt.setMt(x(0));
    AlphaS as(oiMt,aSMZ);

    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMt,pdg2014::Gf, as(oiMt.Mt()), oiMt.Mt(), order::all);
    
    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    
    long double mu2 = pow(mu,2);
    
    lout(logDEBUG) << "Zero of \\lambda at Lambda_I = " << mu;

    // quadratical difference from zero
    state_type a = avP2MS(mu2);
    fvec(0) = a[couplings::lam];

    lout(logDEBUG) << "lam  = " << a[couplings::lam];
    return 0;
  }
};

// Meta-stability - Stability bound
struct Mt_Beta0 : Functor<double>
{
  OSinput oi;
  long double mu;
  long double aSMZ;

  Mt_Beta0(void): Functor<double>(1,1) {}

  // Explicit dependence on asMZ is usefull for plots 
  // depending on asMZ errors
  Mt_Beta0(const OSinput& in_, long double mu_, long double alphaS_ = pdg2014::asMZ) : oi(in_), mu(mu_), aSMZ(alphaS_), Functor<double>(1,1)
  {
  }

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {
    OSinput oiMt(oi);
    oiMt.setMt(x(0));
    
    AlphaS as(oiMt,aSMZ);

    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMt,pdg2014::Gf, as(oiMt.Mt()), oiMt.Mt(), order::all);
    
    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    
    long double mu2 = pow(mu,2);
    
    lout(logDEBUG) << "Zero of \\lambda at Lambda_I = " << mu;

    // quadratical difference from zero
    std::pair<state_type, state_type> ab = avP2MS.AandB(mu2);
    fvec(0) = ab.second[couplings::lam]; // beta = 0

    lout(logDEBUG) << "\\beta_lam  = " << ab.second[couplings::lam];
    return 0;
  }
};

// Meta-stability - Stability bound
struct MH_Beta0 : Functor<double>
{
  OSinput oi;
  long double asMt;
  long double mu;

  MH_Beta0(void): Functor<double>(1,1) {}

  // Explicit dependence on asMZ is usefull for plots 
  // depending on asMZ errors
  MH_Beta0(const OSinput& in_, long double mu_, long double alphaS = pdg2014::asMZ) : oi(in_), mu(mu_), Functor<double>(1,1)
  {
    AlphaS as(oi,alphaS);
    asMt = as(oi.Mt());
  }

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {

    OSinput oiMH(oi);
    // double Mt = x(0);
    oiMH.setMH(x(0));
    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMH,pdg2014::Gf, asMt, oiMH.Mt(), order::all);
    
    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    
    long double mu2 = pow(mu,2);
    
    lout(logDEBUG) << "Zero of \\lambda at Lambda_I = " << mu;

    // quadratical difference from zero
    std::pair<state_type, state_type> ab = avP2MS.AandB(mu2);
    fvec(0) = ab.second[couplings::lam]; // beta = 0

    lout(logDEBUG) << "\\beta_lam  = " << ab.second[couplings::lam];
    return 0;
  }
};


long double critMH(const OSinput&, long double);

long double critMt(const OSinput&, long double);


long double critMH2(const OSinput&, long double);

long double critMH_scaleNotFixed(const OSinput&, long double);

long double critMt_test(const OSinput& oi, long double mu, long double alphaS);

template<typename BoundT>
std::pair<long double, int> critMH0(const OSinput&, long double, long double);

template<typename BoundT>
std::pair<long double, int> critMt0(const OSinput&, long double, long double);

long double critMt(const OSinput&, long double);

template<typename BoundT>
std::pair<long double,long double> critMt_scaleNotFixed(const OSinput& oi, long double mu, long double alphaS = pdg2014::asMZ);

#endif  // __TOOLS_HPP__
