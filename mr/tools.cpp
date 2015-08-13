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
  
  SMCouplings::value_type operator()(long double MH)
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

struct bLambdaAtMuMH
{
  
  OSinput oi;
  long double mu2;
  long double asMt;
  
  bLambdaAtMuMH(OSinput in_, long double mu2_) : oi(in_), mu2(mu2_)
  {
    AlphaS as(oi);
    asMt = as(oi.Mt());
  }
  
  SMCouplings::value_type operator()(MRt MH)
  {
    oi.setMH(MH);
    // Set of all running parameters at scale Mt
    P2MS pMSmt(oi,pdg2014::Gf, asMt, oi.Mt(), order::all);

    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    return avP2MS.AandB(mu2).second[couplings::lam];
  }
  typedef SMCouplings::value_type value_type;
};


MRt critMH(const OSinput& oi, long double mu)
{
  bLambdaAtMuMH lambdaMPL(oi, pow(mu,2));
  tolerance tol = 1e-8;
  
  // Finding bounds with different lambda sign
  SMCouplings::value_type startLam = lambdaMPL(oi.MH());
  MRt leftMH, rightMH;
  // If initial lambda is negative which is true for most of realistic SM
  // inputs than we should increase MH by 10GeV before we finish with
  // positive, but not rising to Landau pole lambda

  // if(startLam < 0)
  //   {
  //     // We already know lower bound on Higgs mass
  //     leftMH = oi.MH();

  //     rightMH = oi.MH() + 10.;  // 10 GeV larger
  //     while(lambdaMPL(rightMH) < 0) rightMH += 10.;
  //   }
  // else
  //   {
  //     // We already know upper bound on Higgs mass
  //     rightMH = oi.MH();
      
  //     leftMH = oi.MH() - 10.;  // 10 GeV larger
  //     while(lambdaMPL(leftMH) > 0) leftMH -= 10.;
  //   }

  leftMH = 120;
  rightMH = 130;
  std::cout <<"Bisection for Higgs mass in interval [" << leftMH << ", " << rightMH << "]" << std::endl;
  std::pair<MRt,MRt> MHcritInterval = boost::math::tools::bisect(lambdaMPL, leftMH, rightMH, tol);
  
  boost::numeric::interval<MRt> MHint(MHcritInterval.first, MHcritInterval.second);
  
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
  }
  
  SMCouplings::value_type operator()(long double Mt)
  {

    oi.setMt(Mt);
    AlphaS as(oi);
    asMt = as(oi.Mt());


    // Set of all running parameters at scale Mt
    P2MS pMSmt(oi,pdg2014::Gf, asMt, oi.Mt(), order::all);

    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    return avP2MS(mu2)[couplings::lam];
  }

};

struct bLambdaAtMuMt
{
  
  OSinput oi;
  long double mu2;
  long double asMt;
  
  bLambdaAtMuMt(OSinput in_, long double mu2_) : oi(in_), mu2(mu2_)
  {
  }
  
  SMCouplings::value_type operator()(long double Mt)
  {
    oi.setMt(Mt);
    AlphaS as(oi);
    asMt = as(oi.Mt());

    // Set of all running parameters at scale Mt
    P2MS pMSmt(oi,pdg2014::Gf, asMt, oi.Mt(), order::all);

    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    return avP2MS.AandB(mu2).second[couplings::lam];
  }

};



long double critMt(const OSinput& oi, long double mu)
{
  bLambdaAtMuMt lambdaMPL(oi, pow(mu,2));
  tolerance tol = 1e-8;

  // Finding bounds with different lambda sign
  SMCouplings::value_type startLam = lambdaMPL(oi.Mt());
  long double leftMt, rightMt;
  // If initial lambda is negative which is true for most of realistic SM
  // inputs than we should decrease Mt by 10GeV before we finish with
  // positive, but not rising to Landau pole lambda

  // For larger Mt beta_lam(\mu)=0 for smaller \mu

  // if(startLam < 0)
  //   {

  //     std::cout << "Starting value for lambda is negative, decreasing Top mass" << std::endl;
  //     // We already know upper bound on Top mass
  //     rightMt = oi.Mt();

  //     leftMt = oi.Mt() - 10.;  // 10 GeV larger
  //     while(lambdaMPL(leftMt) < 0) leftMt -= 10.;
  //   }
  // else
  //   {

  //     std::cout << "Starting value for lambda is positive, increasing Top mass" << std::endl;
  //     // We already know upper bound on Higgs mass
  //     leftMt = oi.Mt();
      
  //     rightMt = oi.Mt() + 10.;  // 10 GeV larger
  //     while(lambdaMPL(rightMt) > 0) rightMt += 10.;
  //   }
  
  leftMt = 169;
  rightMt = 171;

  std::cout <<"\n\n*******\n\tBisection for Top mass in interval [" << leftMt << ", " << rightMt << "]" << std::endl << std::endl << std::endl;
  std::pair<long double, long double> MtcritInterval = boost::math::tools::bisect(lambdaMPL, leftMt, rightMt, tol);
  
  boost::numeric::interval<long double> Mtint(MtcritInterval.first, MtcritInterval.second);
  
  std::cout << std::setprecision(10);
  std::cout << "==> critical Mt  at [lambda(Mpl)=0] = [" << MtcritInterval.first << ',' << MtcritInterval.second << "]\n";
  

  return boost::numeric::median(Mtint);

}



struct MH_functor : Functor<MRt>
{
  OSinput oi;
  long double mu2;
  long double asMt;
  
  MH_functor(void): Functor<MRt>(1,1) {}

  MH_functor(const OSinput& in_, long double mu2_) : oi(in_), mu2(mu2_), Functor<MRt>(1,1)
  {
    AlphaS as(oi);
    asMt = as(oi.Mt());
  }

  int operator()(const EigenMVec &x, EigenCVec &fvec) const
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
  
  
  EigenMVec x(1);
  x(0) = oi.MH();
  std::cout << "x: " << x << std::endl;
  
  MH_functor functor(oi, pow(mu,2));
  Eigen::NumericalDiff<MH_functor> numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MH_functor>, MRt> lm(numDiff);
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




struct MH_functor2 : Functor<MRt>
{
  OSinput oi;
  long double asMt;
  
  MH_functor2(void): Functor<MRt>(2,2) {}

  MH_functor2(const OSinput& in_) : oi(in_),Functor<MRt>(2,2)
  {
    AlphaS as(oi);
    asMt = as(oi.Mt());
  }

  int operator()(const EigenMVec &x, EigenCVec &fvec) const
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
    std::pair<SMCouplings, SMCouplings> ab = avP2MS.AandB(mu2);
    fvec(0) = pow(ab.first[couplings::lam],2);
    fvec(1) = pow(ab.second[couplings::lam],2);

    return 0;
  }
};




long double critMH_scaleNotFixed(const OSinput& oi, long double mu)
{
  
  EigenMVec x(2);
  x(0) = oi.MH();
  x(1) = log10(pow(mu,2));
  std::cout << "x: " << x << std::endl;
  
  MH_functor2 functor(oi);
  Eigen::NumericalDiff<MH_functor2> numDiff(functor);
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<MH_functor2>,MRt> lm(numDiff);
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





// 
//    <BoundT> - is functor type for boundary condition
// 
//    <Mt_Stability is> criterium for absolute stability
//    lambda=beta_lambda=0
// 
//    <Mt_Instability> corresponds to border between instability 
//    and meta-stability regions
//    beta_lambda=0, lambda=lambda_min, where lambda_min is
//   such that tunneling rate P(lambda)=1
// 
// template <typename BoundT>
// std::pair<long double,long double> critMt_scaleNotFixed(const OSinput& oi, long double mu, long double alphaS)
// {
  
  
//   Eigen::VectorXd x(2);
//   x(0) = oi.Mt();
//   x(1) = log10(pow(mu,2));
  
//   BoundT functor(oi, alphaS);
//   Eigen::NumericalDiff<BoundT> numDiff(functor);
//   Eigen::LevenbergMarquardt<Eigen::NumericalDiff<BoundT>,double> lm(numDiff);
//   lm.parameters.maxfev = 200;
//   lm.parameters.xtol = 1.0e-12;
//   lm.parameters.epsfcn = 1.0e-5;

//   int ret = lm.minimize(x);
//   lout(logINFO) << "Number of function evaluations" << lm.iter;
//   lout(logINFO) << "Status of minimization        " << ret;

//   if (ret != 2) 
//     lout(logERROR) << "Minimization problems, status: " << ret;

//   lout(logINFO)  << "Mt           = " << x(0);
//   lout(logINFO)  << "log10(mu^2)  = " << x(1);
//   lout(logDEBUG) << "lambda       = " << lm.fvec(0);
//   lout(logDEBUG) << "beta(lambda) = " << lm.fvec(1);
//   // lout(logDEBUG) << "Mt(f)        = " << lm.fvec(2);
  
//   return std::pair<long double,long double>(x(0), x(1));
// }


// Meta-stability - Stability bound
struct Mt_pointFunctor : Functor<MRt>
{
  OSinput oi;
  long double asMt;
  
  Mt_pointFunctor(void): Functor<MRt>(1,2) {}

  // Explicit dependence on asMZ is usefull for plots 
  // depending on asMZ errors
  Mt_pointFunctor(const OSinput& in_, long double alphaS = pdg2014::asMZ) : oi(in_),Functor<MRt>(1,2)
  {
    AlphaS as(oi,alphaS);
    asMt = as(oi.Mt());
  }

  int operator()(const EigenMVec &x, EigenCVec &fvec) const
  {


    OSinput oiMt(oi);
    // double Mt = x(0);
    // oiMt.setMt(x(0));
    // Set of all running parameters at scale Mt
    P2MS pMSmt(oiMt,pdg2014::Gf, asMt, oiMt.Mt(), order::all);
    
    Couplings<3,3,3,
              3,3,-1,
              3,0,0> avP2MS(pMSmt);

    
    long double mu2 = pow(10.,x(0));
    
    lout(logDEBUG) << "Instability scale Lambda_I = " << sqrt(mu2);

    // quadratical difference from zero
    std::pair<SMCouplings, SMCouplings> ab = avP2MS.AandB(mu2);
    fvec(0) = pow(ab.first[couplings::lam],1);
    fvec(1) = pow(ab.second[couplings::lam],1);

    lout(logDEBUG) << "lam  = " << ab.first[couplings::lam];
    lout(logDEBUG) << "blam = " << ab.second[couplings::lam];
    // We are interested in minimal Mt 
    // corresponding to border of Stability
    // fvec(2) = pow(x(0),2);
    return 0;
  }
};


long double critMt_test(const OSinput& oi, long double mu, long double alphaS)
{
  
  
  EigenMVec x(1);
  x(0) = log10(pow(mu,2));
  
  Mt_pointFunctor functor(oi, alphaS);

  Eigen::HybridNonLinearSolver<Mt_pointFunctor, MRt> solver(functor);
  int ret = solver.hybrd1(x);

  // Eigen::NumericalDiff<BoundT> numDiff(functor);
  // Eigen::LevenbergMarquardt<Eigen::NumericalDiff<BoundT>,double> lm(numDiff);
  // solver.parameters.maxfev = 200;
  // solver.parameters.xtol = 1.0e-12;
  // solver.parameters.epsfcn = 1.0e-5;

  // int ret = lm.minimize(x);
  lout(logINFO) << "Number of function evaluations" << solver.iter;
  lout(logINFO) << "Status of minimization        " << ret;

  if (ret != 2) 
    lout(logERROR) << "Minimization problems, status: " << ret;

  lout(logINFO)  << "log10(mu^2)  = " << x(0);
  lout(logDEBUG) << "lambda       = " << solver.fvec(0);
  // lout(logDEBUG) << "beta(lambda) = " << solver.fvec(1);
  // lout(logDEBUG) << "Mt(f)        = " << lm.fvec(2);
  
  return x(0);
}


template <typename BoundT>
std::pair<long double,long double> critMt_scaleNotFixed(const OSinput& oi, long double mu, long double alphaS)
{
  
  
  EigenMVec x(2);
  x(0) = oi.Mt();
  x(1) = log10(pow(mu,2));
  
  BoundT functor(oi, alphaS);

  Eigen::HybridNonLinearSolver<BoundT, MRt> solver(functor);
  int ret = solver.hybrd1(x);

  lout(logINFO) << "Number of function evaluations" << solver.iter;
  lout(logINFO) << "Status of minimization        " << ret;

  if (ret != 2) 
    lout(logERROR) << "Minimization problems, status: " << ret;

  lout(logINFO)  << "Mt           = " << x(0);
  lout(logINFO)  << "log10(mu^2)  = " << x(1);
  lout(logDEBUG) << "lambda       = " << solver.fvec(0);
  lout(logDEBUG) << "beta(lambda) = " << solver.fvec(1);
  
  return std::pair<long double,long double>(x(0), x(1));
}

template std::pair<long double,long double> critMt_scaleNotFixed<Mt_Stability>(const OSinput&, long double, long double);
// template std::pair<long double,long double> critMt_scaleNotFixed<Mt_StabilityScale>(const OSinput&, long double, long double);
template std::pair<long double,long double> critMt_scaleNotFixed<Mt_Instability>(const OSinput&, long double, long double);

template<typename BoundT>
std::pair<long double, int> critMH0(const OSinput& oi, long double mu, long double aSMZ)
{
  
  EigenMVec x(1);
  x(0) = oi.MH();

  BoundT functor(oi, mu, aSMZ);

  Eigen::HybridNonLinearSolver<BoundT, MRt> solver(functor);
  solver.diag.setConstant(1, 1.);
  solver.useExternalScaling = true;
  solver.parameters.xtol = 1.0e-12;
  solver.parameters.epsfcn = 1.0e-5;

  int ret = solver.solveNumericalDiff(x);

  lout(logINFO) << "Number of function evaluations" << solver.iter;
  lout(logINFO) << "Status of minimization        " << ret;
  lout(logINFO) << "MH           = " << x(0);
  
  if (ret != 2) 
    lout(logERROR) << "Minimization problems, status: " << ret;
  
  return std::make_pair(x(0), ret);
}


template<typename BoundT>
std::pair<long double, int> critMt0(const OSinput& oi, long double mu, long double aSMZ)
{
  
  EigenMVec x(1);
  x(0) = oi.Mt();

  BoundT functor(oi, mu, aSMZ);

  Eigen::HybridNonLinearSolver<BoundT, MRt> solver(functor);
  solver.diag.setConstant(1, 1.);
  solver.useExternalScaling = true;
  solver.parameters.xtol = 1.0e-12;
  solver.parameters.epsfcn = 1.0e-5;

  int ret = solver.solveNumericalDiff(x);

  lout(logINFO) << "Number of function evaluations" << solver.iter;
  lout(logINFO) << "Status of minimization        " << ret;
  lout(logINFO) << "Mt           = " << x(0);
  
  if (ret != 2) 
    lout(logERROR) << "Minimization problems, status: " << ret;
  
  return std::make_pair(x(0), ret);
}

// H
// template std::pair<long double, int> critMH0<MH_Lambda0>(const OSinput&, long double, long double);
template std::pair<long double, int> critMH0<MH_Beta0>(const OSinput&, long double, long double);

// t
template std::pair<long double, int> critMt0<Mt_Lambda0>(const OSinput&, long double, long double );
template std::pair<long double, int> critMt0<Mt_Beta0>(const OSinput&, long double, long double);

loglevel_e loglevel = logERROR;
