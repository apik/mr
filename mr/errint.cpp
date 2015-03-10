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

#include "errint.hpp"
#include "base.hpp"
#include "WW.hpp"
#include <bitset>
#include <iostream>
#include <list>
WWerr::WWerr(OSinputErr oierr, long double mu)
{

  std::list<PoleMass*> wl;

  for(size_t mbi=0; mbi < 2; mbi++)
      for(size_t mWi=0; mWi < 2; mWi++)
          for(size_t mZi=0; mZi < 2; mZi++)
              for(size_t mHi=0; mHi < 2; mHi++)
                  for(size_t mti=0; mti < 2; mti++)
                    {
                      OSinput oi(
                                 mbi == 0? oierr.Mb().lower():oierr.Mb().upper(),
                                 mWi == 0? oierr.MW().lower():oierr.MW().upper(),
                                 mZi == 0? oierr.MZ().lower():oierr.MZ().upper(),
                                 mHi == 0? oierr.MH().lower():oierr.MH().upper(),
                                 mti == 0? oierr.Mt().lower():oierr.Mt().upper()
                                 );
                      // WW<OS>* wwtmp = new WW<OS>(oi, 173*173)
                        wl.push_back(new WW<OS>(oi, 173*173));
                    }
}
