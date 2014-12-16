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

#include "mfitter.hpp"
#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp"
#include "dr.hpp"

std::vector<long double> observables(long double gp,
                                     long double g,
                                     long double gs,
                                     long double yt,
                                     long double lam,
                                     long double mu0,
                                     long double mu2)
{

  MSinput mi = MSinput::fromConsts(mu2, mu0, lam, 0, yt, g, gp);

  std::cout << "In fitter:" << std::endl;
  std::cout << mi.mW() << std::endl;
  std::cout << mi.mZ() << std::endl;
  std::cout << mi.mH() << std::endl;
  std::cout << mi.mt() << std::endl;
  std::cout << mi.vev() << std::endl;
  
  WW<MS> MWth(mi, mu2);
  ZZ<MS> MZth(mi, mu2);
  HH<MS> MHth(mi, mu2);
  tt<MS> Mtth(mi, mu2);
  dr<MS> DRth(mi, mu2);

  long double aEW  = mi.alpha()/4./Pi;
  long double aQCD = pow(gs/4./Pi,2);
  
  long double MWpole = mi.mW()*sqrt(1 +
                                    aEW*MWth.x10() +
                                    aEW*aQCD*MWth.x11() +
                                    aEW*aEW*MWth.x20());
  
  long double MZpole = mi.mZ()*sqrt(1 +
                                    aEW*MZth.x10() +
                                    aEW*aQCD*MZth.x11() +
                                    aEW*aEW*MZth.x20());
  
  long double MHpole = mi.mH()*sqrt(1 +
                                    aEW*MHth.x10() +
                                    aEW*aQCD*MHth.x11() +
                                    aEW*aEW*MHth.x20());
  
  long double Mtpole = mi.mt()*(1 +
                                aEW*Mtth.x10() +
                                aEW*aQCD*Mtth.x11() +
                                aEW*aEW*Mtth.x20() +
                                // QCD
                                aQCD*Mtth.x01() +
                                aQCD*aQCD*Mtth.x02() +
                                aQCD*aQCD*aQCD*Mtth.x03()
                                );
  
  long double Gf = (1 +
                    aEW*DRth.dr10() +
                    aEW*aQCD*DRth.dr11() +
                    aEW*aEW*DRth.dr20())/sqrt(2)/pow(mi.vev(),2);
  
  std::vector<long double> obs(5);
  
  obs[0] = MWpole;
  obs[1] = MZpole;
  obs[2] = MHpole;
  obs[3] = Mtpole;
  obs[4] = Gf;

  return obs;
}


double MfitterFcn::operator()(const std::vector<double>& par) const {
  
  // 6 parameters: gp, g, gs, yt, lambda, mu0
  // assert(par.size() == 6);

  std::vector<long double> ob = observables(par[0], // 0.35761
                                            par[1], // 0.64822, 
                                            par[2], // 1.1666,
                                            par[3], // 0.93558,
                                            par[4], // 0.12711,
                                            par[5], // 132.03,
                                            pow(theScale,2)
                                            );

  ob.push_back(par[2]*par[2]/4./Pi);               // alpha_S =g3^2/4/Pi

  // GaussFunction gauss(par[0], par[1], par[2]);

  long double chi2 = 0.;
  for(unsigned int n = 0; n < theMeasurements.size(); n++) {
    chi2 += ((ob[n] - theMeasurements[n])*(ob[n] - theMeasurements[n])/theMVariances[n]);
    std::cout << "obs = " << ob[n] << "  and mes = " << theMeasurements[n] << std::endl;
  }

  return chi2;
}
