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

#include "boost/math/tools/roots.hpp"
#include "boost/numeric/interval.hpp"
#include "tools.hpp"
#include "p2ms.hpp"

#include "betaQCD.hpp"
#include "betaSM.hpp"



class tolerance {
public:
  tolerance(long double eps) :
    _eps(eps) {
  }
  bool operator()(long double a, long double b) {
    return (fabs(b - a) <= _eps);
  }
private:
  long double _eps;
};



struct LambdaAtMuMH
{
  
  OSinput oi;
  long double mu2;
  long double asMt;
  
  LambdaAtMuMH(OSinput in_, long double mu2_) : oi(in_), mu2(mu2_)
  {
    AlphaS as(oi);
    asMt = as(oi.Mt());
  }
  
  long double operator()(long double MH)
  {
    oi.setMH(MH);
    // Set of all running parameters at scale Mt
    P2MS pMSmt(oi,pdg2014::Gf, asMt, oi.Mt(), order::all);

    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    return avP2MS(mu2)[couplings::lam];
  }

};



long double critMH(const OSinput& oi, long double mu)
{
  LambdaAtMuMH lambdaMPL(oi, pow(mu,2));
  tolerance tol = 1e-12;

  // Finding bounds with different lambda sign
  long double startLam = lambdaMPL(oi.MH());
  long double leftMH, rightMH;
  // If initial lambda is negative which is true for most of realistic SM
  // inputs than we should increase MH by 10GeV before we finish with
  // positive, but not rising to Landau pole lambda

  if(startLam < 0)
    {
      // We already know lower bound on Higgs mass
      leftMH = oi.MH();

      rightMH = oi.MH() + 10.;  // 10 GeV larger
      while(lambdaMPL(rightMH) < 0) rightMH += 10.;
    }
  else
    {
      // We already know upper bound on Higgs mass
      rightMH = oi.MH();
      
      leftMH = oi.MH() - 10.;  // 10 GeV larger
      while(lambdaMPL(leftMH) > 0) leftMH -= 10.;
    }

  std::cout <<"Bisection for Higgs mass in interval [" << leftMH << ", " << rightMH << "]" << std::endl;
  std::pair<long double, long double> MHcritInterval = boost::math::tools::bisect(lambdaMPL, leftMH, rightMH, tol);
  
  boost::numeric::interval<long double> MHint(MHcritInterval.first, MHcritInterval.second);
  
  std::cout << std::setprecision(10);
  std::cout << "==> critical MH  at [lambda(Mpl)=0] = [" << MHcritInterval.first << ',' << MHcritInterval.second << "]\n";
  

  return boost::numeric::median(MHint);

}

struct LambdaAtMuMt
{
  
  OSinput oi;
  long double mu2;
  long double asMt;
  
  LambdaAtMuMt(OSinput in_, long double mu2_) : oi(in_), mu2(mu2_)
  {
    AlphaS as(oi);
    asMt = as(oi.Mt());
  }
  
  long double operator()(long double Mt)
  {
    oi.setMt(Mt);
    // Set of all running parameters at scale Mt
    P2MS pMSmt(oi,pdg2014::Gf, asMt, oi.Mt(), order::all);

    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    return avP2MS(mu2)[couplings::lam];
  }

};



long double critMt(const OSinput& oi, long double mu)
{
  LambdaAtMuMt lambdaMPL(oi, pow(mu,2));
  tolerance tol = 1e-12;

  // Finding bounds with different lambda sign
  long double startLam = lambdaMPL(oi.Mt());
  long double leftMt, rightMt;
  // If initial lambda is negative which is true for most of realistic SM
  // inputs than we should decrease Mt by 10GeV before we finish with
  // positive, but not rising to Landau pole lambda

  // For larger Mt beta_lam(\mu)=0 for smaller \mu
  if(startLam < 0)
    {

      std::cout << "Starting value for lambda is negative, decreasing Top mass" << std::endl;
      // We already know upper bound on Top mass
      rightMt = oi.Mt();

      leftMt = oi.Mt() - 10.;  // 10 GeV larger
      while(lambdaMPL(leftMt) < 0) leftMt -= 10.;
    }
  else
    {

      std::cout << "Starting value for lambda is positive, increasing Top mass" << std::endl;
      // We already know upper bound on Higgs mass
      leftMt = oi.Mt();
      
      rightMt = oi.Mt() + 10.;  // 10 GeV larger
      while(lambdaMPL(rightMt) > 0) rightMt += 10.;
    }

  std::cout <<"\n\n*******\n\tBisection for Top mass in interval [" << leftMt << ", " << rightMt << "]" << std::endl << std::endl << std::endl;
  std::pair<long double, long double> MtcritInterval = boost::math::tools::bisect(lambdaMPL, leftMt, rightMt, tol);
  
  boost::numeric::interval<long double> Mtint(MtcritInterval.first, MtcritInterval.second);
  
  std::cout << std::setprecision(10);
  std::cout << "==> critical Mt  at [lambda(Mpl)=0] = [" << MtcritInterval.first << ',' << MtcritInterval.second << "]\n";
  

  return boost::numeric::median(Mtint);

}



struct MH_functor : Functor<double>
{
  OSinput oi;
  long double mu2;
  long double asMt;
  
  MH_functor(void): Functor<double>(1,1) {}

  MH_functor(const OSinput& in_, long double mu2_) : oi(in_), mu2(mu2_), Functor<double>(1,1)
  {
    AlphaS as(oi);
    asMt = as(oi.Mt());
  }

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {


    OSinput oiMH(oi);
    // double MH = x(0);
    oiMH.setMH(x(0));
    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMH,pdg2014::Gf, asMt, oiMH.Mt(), order::all);
    
    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    // quadratical difference from zero
    fvec(0) = pow(avP2MS(mu2)[couplings::lam],2);

    return 0;
  }
};




long double critMH2(const OSinput& oi, long double mu)
{
  
  
  Eigen::VectorXd x(1);
  x(0) = oi.MH();
  std::cout << "x: " << x << std::endl;
  
  MH_functor functor(oi, pow(mu,2));
  Eigen::NumericalDiff<MH_functor> numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MH_functor>,double> lm(numDiff);
  lm.parameters.maxfev = 200;
  lm.parameters.xtol = 1.0e-12;
  lm.parameters.epsfcn = 1.0e-5;
  std::cout << "Max fev= " << lm.parameters.maxfev << std::endl;

  int ret = lm.minimize(x);
  std::cout << lm.iter << std::endl;
  std::cout << "Minimum " << ret << std::endl;
  
  std::cout << "x that minimizes the function: " << x << std::endl;

  return x(0);

}




struct MH_functor2 : Functor<double>
{
  OSinput oi;
  long double asMt;
  
  MH_functor2(void): Functor<double>(2,2) {}

  MH_functor2(const OSinput& in_) : oi(in_),Functor<double>(2,2)
  {
    AlphaS as(oi);
    asMt = as(oi.Mt());
  }

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {


    OSinput oiMH(oi);
    // double MH = x(0);
    oiMH.setMH(x(0));
    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMH,pdg2014::Gf, asMt, oiMH.Mt(), order::all);
    
    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    
    long double mu2 = pow(10.,x(1));
    
    // quadratical difference from zero
    std::pair<state_type, state_type> ab = avP2MS.AandB(mu2);
    fvec(0) = pow(ab.first[couplings::lam],2);
    fvec(1) = pow(ab.second[couplings::lam],2);

    return 0;
  }
};




long double critMH_scaleNotFixed(const OSinput& oi, long double mu)
{
  
  Eigen::VectorXd x(2);
  x(0) = oi.MH();
  x(1) = log10(pow(mu,2));
  std::cout << "x: " << x << std::endl;
  
  MH_functor2 functor(oi);
  Eigen::NumericalDiff<MH_functor2> numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MH_functor2>,double> lm(numDiff);
  lm.parameters.maxfev = 200;
  lm.parameters.xtol = 1.0e-13;
  lm.parameters.epsfcn = 1.0e-6;
  std::cout << "Max fev= " << lm.parameters.maxfev << std::endl;

  int ret = lm.minimize(x);
  std::cout << "Affter " << lm.iter << " iterations," << std::endl;
  std::cout << "Status: " << ret << std::endl;
  
  std::cout << "x that minimizes the function: " << x << std::endl;

  return x(0);

}









struct Mt_functor2 : Functor<double>
{
  OSinput oi;
  long double asMt;
  
  Mt_functor2(void): Functor<double>(2,2) {}

  Mt_functor2(const OSinput& in_) : oi(in_),Functor<double>(2,2)
  {
    AlphaS as(oi);
    asMt = as(oi.Mt());
  }

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {


    OSinput oiMt(oi);
    // double MH = x(0);
    oiMt.setMt(x(0));
    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMt,pdg2014::Gf, asMt, oiMt.Mt(), order::all);
    
    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    
    long double mu2 = pow(10.,x(1));
    
    // quadratical difference from zero
    std::pair<state_type, state_type> ab = avP2MS.AandB(mu2);
    fvec(0) = pow(ab.first[couplings::lam],2);
    fvec(1) = pow(ab.second[couplings::lam],2);

    return 0;
  }
};




long double critMt_scaleNotFixed(const OSinput& oi, long double mu)
{
  
  
  Eigen::VectorXd x(2);
  x(0) = oi.Mt();
  x(1) = log10(pow(mu,2));
  std::cout << "x: " << x << std::endl;
  
  Mt_functor2 functor(oi);
  Eigen::NumericalDiff<Mt_functor2> numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<Mt_functor2>,double> lm(numDiff);
  lm.parameters.maxfev = 200;
  lm.parameters.xtol = 1.0e-13;
  lm.parameters.epsfcn = 1.0e-6;
  std::cout << "Max fev= " << lm.parameters.maxfev << std::endl;

  int ret = lm.minimize(x);
  std::cout << lm.iter << std::endl;
  std::cout << "Minimum " << ret << std::endl;
  
  std::cout << "x that minimizes the function: " << x << std::endl;

  return x(0);

}

