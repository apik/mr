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

#include <bb.hpp>
#include <tt.hpp>
//#include "CRunDec.h"

long double mOS2mMS(long double MM, long double mu2, size_t nl, size_t nh, size_t loops = 1)
{
  double NL = double(nl), NH = double (nh);
  if(loops == 1)
    {
      long double mOS2mMS1 =   + log(mu2/MM) * (  - 4 );
      
      mOS2mMS1 +=  - 16./3.;

      return mOS2mMS1;
    }
  else if(loops == 2)
    {
      long double mOS2mMS2 =   + log(mu2/MM) * (  - 314./3. + 52./9.*NL + 52./9.*NH );

      mOS2mMS2 +=  + pow(log(mu2/MM),2) * (  - 14 + 4./3.*NL + 4./3.*NH );

      mOS2mMS2 +=  + Zeta2 * (  - 64./3. - 32./3.*log(2.) + 16./3.*NL - 32./3.*NH );

      mOS2mMS2 +=  - 3305./18. + 8./3.*Zeta3 + 71./9.*NL + 143./9.*NH;

      return mOS2mMS2;
    }
  else if(loops == 3)
    {
      long double mOS2mMS3 = + log(mu2/MM) * (  - 42650./9. + 48*Zeta3 + 14164./27.*NL + 448./9.*NL
                          *Zeta3 - 712./81.*pow(NL,2) + 18052./27.*NH + 448./9.*NH*Zeta3 - 
                          2288./81.*NH*NL - 1576./81.*pow(NH,2) );
      
      mOS2mMS3 +=  + pow(log(mu2/MM),2) * (  - 3034./3. + 428./3.*NL - 104./27.*pow(NL,2) 
                                   + 428./3.*NH - 208./27.*NH*NL - 104./27.*pow(NH,2) );
      
      mOS2mMS3 +=  + pow(log(mu2/MM),3) * (  - 84 + 128./9.*NL - 16./27.*pow(NL,2) + 
                                   128./9.*NH - 32./27.*NH*NL - 16./27.*pow(NH,2) );
      
      mOS2mMS3 +=  + Zeta2 * (  - 99980./27. - 4928./3.*log(2.) + 896./9.*pow(log(2.),2) 
                              + 10648./9.*Zeta3 + 15056./27.*NL + 1408./27.*NL*log(2.) - 
                              256./27.*NL*pow(log(2.),2) - 416./27.*pow(NL,2) - 215728./81.*NH
                              + 81920./27.*NH*log(2.) + 128./27.*NH*pow(log(2.),2) + 96*NH*Zeta3 + 
                              416./27.*NH*NL + 512./135.*pow(NH,2) );

      mOS2mMS3 +=  + Zeta2*log(mu2/MM) * (  - 384 - 192*log(2.) + 1120./9.*NL + 128./9.*NL*
                                log(2.) - 64./9.*pow(NL,2) - 1472./9.*NH + 128./9.*NH*log(2.) + 64./9.
                                *NH*NL + 128./9.*pow(NH,2) );

      mOS2mMS3 +=  - 1259285./162. + 4864./9.*a4 + 608./27.*pow(log(2.),4) - 
        13640./27.*Zeta5 + 6820./9.*Zeta4 + 584./9.*Zeta3 + 172318./243.*NL
          - 512./27.*NL*a4 - 64./81.*NL*pow(log(2.),4) - 4880./27.*NL*Zeta4
        + 5656./27.*NL*Zeta3 - 4706./729.*pow(NL,2) - 224./27.*pow(NL,2)*Zeta3 
        + 315526./243.*NH - 512./27.*NH*a4 - 64./81.*NH*pow(log(2.),4) 
        - 80*NH*Zeta5 - 6560./27.*NH*Zeta4 - 6008./27.*NH*Zeta3 - 
        23668./729.*NH*NL + 128./27.*NH*NL*Zeta3 - 18962./729.*pow(NH,2)
        + 352./27.*pow(NH,2)*Zeta3;
      
      return mOS2mMS3;
    }
  else 
    return 0;
}

std::complex<long double> bb::x02(size_t nL, size_t nH, size_t boson)
{     
  return std::complex<long double>(mOS2mMS(MMb, mu2, nL, nH, 2), 0);
}

std::complex<long double> bb::x03(size_t nL, size_t nH, size_t boson)
{     
  return std::complex<long double>(mOS2mMS(MMb, mu2, nL, nH, 3), 0);
}

std::complex<long double> tt::x02(size_t nL, size_t nH, size_t boson)
{     
  return std::complex<long double>(mOS2mMS(MMb, mu2, nL + nH, nH, 2), 0);
}

std::complex<long double> tt::x03(size_t nL, size_t nH, size_t boson)
{     
  return std::complex<long double>(mOS2mMS(MMb, mu2, nL + nH, nH, 3), 0);
}
