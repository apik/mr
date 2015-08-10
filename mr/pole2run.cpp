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

#include "pole2run.hpp"
#include <stdexcept>
#include "betaQCD.hpp"
#include "betaQEDQCD.hpp"
#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp" 


RunUpto::RunUpto(OSinput oi, long double al_, long double as_, long double mu_): al0(al_), as0(as_), mu0(mu_)
{
  if(mu0 <  oi.Mb())
    throw std::logic_error("ERROR: input at scale mu > Mb needed");

  // We use matching at mu=Mt
  // And need to run QCD and EW
  // couplings up to Mt first

  AlphaS       as;
  AlphaQEDQCD  al;

  // matching scale at top mass, when full SM aplicable
  ms = oi.MMt();

  WW<OS> w(oi, ms);
  ZZ<OS> z(oi, ms);
  HH<OS> h(oi, ms);
  tt<OS> t(oi, ms);

  long double aQCD = as(ms)/4./Pi;
  long double aEW  = al.QED(ms)/4./Pi;

  // 2-loop EW corrections
  long double dWplus1 = 1 + aEW*w.y10() + aEW*aQCD*w.y11() + aEW*aEW*w.y20();
  long double dZplus1 = 1 + aEW*z.y10() + aEW*aQCD*z.y11() + aEW*aEW*z.y20();
  long double dHplus1 = 1 + aEW*h.y10() + aEW*aQCD*h.y11() + aEW*aEW*h.y20();
  long double dtplus1 = 1 + aEW*t.y10() + aEW*aQCD*t.y11() + aEW*aEW*t.y20()
    + aQCD*t.y01()+ aQCD*aQCD*t.y02()+ aQCD*aQCD*aQCD*t.y03();

  long double Gf = pdg2014::Gf;

  long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
  long double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
    
  a1 = 5./3.*(gg_ggp - gg)/16./Pi/Pi;
  a2 = gg/16./Pi/Pi;
  aS = aQCD;
  ayt = pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
  alam = Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;


  std::cout << " At matching scale mu = " << ms << std::endl;
  std::cout << " g1 = " << sqrt(3./5.*a1)*4*Pi << std::endl;
  std::cout << " g2 = " << sqrt(a2)*4*Pi << std::endl;
  std::cout << " g3 = " << sqrt(aS)*4*Pi << std::endl;
  std::cout << " yt = " << sqrt(ayt)*4*Pi << std::endl;
  std::cout << " lam = " << alam*16*Pi*Pi << std::endl;

  av = new CouplingsSM<3,3,3,3,-1,-1,3>(a1,a2,aS,ayt,0,0,alam,pow(ms,2),3);
  // std::cout << " LAMMMMMMMMMMMM " << av->operator()(pow(10000,2))[6] << " ::: " << lambda(10000) << std::endl;
}

SMCouplings RunUpto::operator()(SMCouplings::value_type mu)
{  
  return av->operator()(pow(mu,2));
}


SMCouplings::value_type RunUpto::lambda(SMCouplings::value_type mu)
{
  return av->operator()(pow(mu,2))[6];
}
