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
#include "p2ms.hpp"
#include <stdexcept>
#include "alphas.hpp"
#include "bb.hpp" 
#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp" 
#include "alphaGF.hpp"

namespace mr
{
  struct DiffGF
  {
    
    OSinput oi;

    unsigned ord;
    double Gf0;
    double alphaS;

    // Precalculated set of corrections
    double Wy10, Wy11, Wy20, Zy10, Zy11, Zy20;

    DiffGF(OSinput in_, double Gf0_, double as_, double mu2, unsigned order_) : oi(in_), Gf0(Gf0_), alphaS(as_), ord(order_)
    {
      WW<OS> dW(oi, mu2);
      ZZ<OS> dZ(oi, mu2);

      Wy10 = dW.y10();
      Wy11 = dW.y11();
      Wy20 = dW.y20();

      Zy10 = dZ.y10();
      Zy11 = dZ.y11();
      Zy20 = dZ.y20();
    }

    double operator()(double alpha)
    {
      double dMyW = 1;
      double dMyZ = 1;
      double Gf;

      if(ord & order::x10)
        {
          dMyW += alpha/4./Pi*Wy10;
          dMyZ += alpha/4./Pi*Zy10;
        }
      if(ord & order::x11)
        {
          dMyW += alpha/4./Pi*alphaS/4./Pi*Wy11;
          dMyZ += alpha/4./Pi*alphaS/4./Pi*Zy11;
        }
      if(ord & order::x20)
        {
          dMyW += pow(alpha/4./Pi,2)*Wy20;
          dMyZ += pow(alpha/4./Pi,2)*Zy20;
        }
    

      Gf = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
      return (Gf0 - Gf)*pow(10,2);
    }

  };
  
  

  
  double AlphaSolve::operator()(const double& mu2 )
  {
    DiffGF dGF(oi, Gf0, alphaS, mu2, ord);

    boost::uintmax_t max_iter=50;
    std::pair<double, double> found = 
      boost::math::tools::toms748_solve(dGF, 1./140., 1./120., tol, max_iter);
  
    boost::numeric::interval<double> fint(found.first, found.second);
  
    std::cout << std::setprecision(10);
    lout(logDEBUG) << "==> 1/alpha = [" << 1./found.first << ',' << 1./found.second << "]";
    
    return boost::numeric::median(fint);;
  }

  
  double AlphaGF::operator()(const double& mu2 )
  {
    
    alphaGF aGF  = alphaGF(oi, mu2);
    
    double daGF = 1;
    
    // Tree level
    double alF = sqrt(2)*Gf0*oi.MMW()/Pi*(1-oi.MMW()/oi.MMZ());
    
    if(ord & order::x10)
      {
        daGF += alF/4./Pi*aGF.a10();
      }
    if(ord & order::x11)
      {
        daGF += alF/4./Pi*alphaS/4./Pi*aGF.a11();
      }
    if(ord & order::x20)
      {
        daGF += pow(alF/4./Pi,2)*aGF.a20();
      }

    lout(logDEBUG) << " Alpha from GF: " << alF*daGF << std::endl;
    return alF*daGF;
  }
  
  // logger
  loglevel_e loglevel = logERROR;

} // namespace mr


