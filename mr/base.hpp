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

#ifndef __BASE_HPP__
#define __BASE_HPP__
#include "logger.hpp"
#include <memory>

namespace mr
{
  class BaseMass
  {
    // 
    // Pure QCD part m_ij=mY_ij by definition
    // 

    // 
    // Mass corrections
    // 
  public:

    virtual long double x10(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;
  
    virtual long double x11(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

    virtual long double x20(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

    virtual long double x01(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Order a^0*as^1 is not implemented for this particle" << std::endl;
      return 0;
    }
    virtual long double x02(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Order a^0*as^2 is not implemented for this particle" << std::endl;
      return 0;
    }
    virtual long double x03(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Order a^0*as^3 is not implemented for this particle" << std::endl;
      return 0;
    }
    virtual long double x04(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Order a^0*as^4 is not implemented for this particle" << std::endl;
      return 0;
    }

    // Gauge-less limit
    virtual long double xgl10(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Gauge-less limit at order a^1*as^0 is not implemented for this particle" << std::endl;
      return 0;
    }
    virtual long double xgl11(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Gauge-less limit at order a^1*as^1 is not implemented for this particle" << std::endl;
      return 0;
    }
    virtual long double xgl20(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Gauge-less limit at order a^2*as^0 is not implemented for this particle" << std::endl;
      return 0;
    }


    // And meta method
    long double x(size_t apow, size_t aspow, size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      if(apow == 1 && aspow == 0)
        return x10(nL, nH, boson);
      if(apow == 1 && aspow == 1)
        return x11(nL, nH, boson);
      if(apow == 2 && aspow == 0)
        return x20(nL, nH, boson);
      if(apow == 0 && aspow == 1)
        return x01(nL, nH, boson);
      if(apow == 0 && aspow == 2)
        return x02(nL, nH, boson);
      if(apow == 0 && aspow == 3)
        return x03(nL, nH, boson);
      if(apow == 0 && aspow == 4)
        return x04(nL, nH, boson);
      return 0;
    }

    long double xgl(size_t apow, size_t aspow, size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      if(apow == 1 && aspow == 0)
        return xgl10(nL, nH, boson);
      if(apow == 1 && aspow == 1)
        return xgl11(nL, nH, boson);
      if(apow == 2 && aspow == 0)
        return xgl20(nL, nH, boson);
      if(apow == 0 && aspow == 1)
        return x01(nL, nH, boson);
      if(apow == 0 && aspow == 2)
        return x02(nL, nH, boson);
      if(apow == 0 && aspow == 3)
        return x03(nL, nH, boson);
      if(apow == 0 && aspow == 4)
        return x04(nL, nH, boson);
      return 0;
    }

  };

  class PoleMassAndCouplings : public BaseMass
  {
  public:
    virtual long double y10(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

    virtual long double y11(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

    virtual long double y20(size_t nL = 2, size_t nH = 1, size_t boson = 1) = 0;

    // Gauge-less limit
    virtual long double ygl10(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Gauge-less limit at order a^1*as^0 is not implemented for this particle" << std::endl;
      return 0;
    }
    virtual long double ygl11(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Gauge-less limit at order a^1*as^1 is not implemented for this particle" << std::endl;
      return 0;
    }
    virtual long double ygl20(size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      lout(logERROR) << "Gauge-less limit at order a^2*as^0 is not implemented for this particle" << std::endl;
      return 0;
    }

    long double y(size_t apow, size_t aspow, size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      if(apow == 1 && aspow == 0)
        return y10(nL, nH, boson);
      if(apow == 1 && aspow == 1)
        return y11(nL, nH, boson);
      if(apow == 2 && aspow == 0)
        return y20(nL, nH, boson);
      if(apow == 0 && aspow == 1)
        return x01(nL, nH, boson);
      if(apow == 0 && aspow == 2)
        return x02(nL, nH, boson);
      if(apow == 0 && aspow == 3)
        return x03(nL, nH, boson);
      if(apow == 0 && aspow == 4)
        return x04(nL, nH, boson);
      return 0;
    }
    
    long double ygl(size_t apow, size_t aspow, size_t nL = 2, size_t nH = 1, size_t boson = 1)
    {
      if(apow == 1 && aspow == 0)
        return ygl10(nL, nH, boson);
      if(apow == 1 && aspow == 1)
        return ygl11(nL, nH, boson);
      if(apow == 2 && aspow == 0)
        return ygl20(nL, nH, boson);
      if(apow == 0 && aspow == 1)
        return x01(nL, nH, boson);
      if(apow == 0 && aspow == 2)
        return x02(nL, nH, boson);
      if(apow == 0 && aspow == 3)
        return x03(nL, nH, boson);
      if(apow == 0 && aspow == 4)
        return x04(nL, nH, boson);
      return 0;
    }


  };
} // namespace mr

#endif  // __BASE_HPP__
