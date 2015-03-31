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
// #include <stdexcept>
#include "betaQCD.hpp"
#include "betaSM.hpp"
// #include "betaQEDQCD.hpp"
// #include "bb.hpp" 
// #include "WW.hpp"
// #include "ZZ.hpp"
// #include "HH.hpp"
// #include "tt.hpp" 



struct LambdaAtMu
{
  
  OSinput oi;
  long double mu2;
  long double asMt;
  
  LambdaAtMu(OSinput in_, long double mu2_) : oi(in_), mu2(mu2_)
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


long double critMH(const OSinput& oi, long double mu)
{
  LambdaAtMu lambdaMPL(oi, pow(mu,2));
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

