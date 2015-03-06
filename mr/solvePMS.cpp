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

#include "solvePMS.hpp"
#include <stdexcept>
#include "betaQCD.hpp"
#include "betaQEDQCD.hpp"
#include "bb.hpp" 
#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp" 


// root-finder
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>



// Difference between calculated Pole masses and input
int Mdiff (const gsl_vector * x, void *params, gsl_vector * f)
{

  OSinput* poi =  ((OSinput *) params);

  // Initial values to fit to
  long double Mb0 = poi->Mb();
  long double MW0 = poi->MW();
  long double MZ0 = poi->MZ();
  long double MH0 = poi->MH();
  long double Mt0 = poi->Mt();
    
  const double g1  = gsl_vector_get (x, 0);
  const double g2  = gsl_vector_get (x, 1);
  const double gs  = 1.1666;
  const double yt  = gsl_vector_get (x, 2);
  const double yb  = gsl_vector_get (x, 3);
  const double lam = gsl_vector_get (x, 4);
  const double muh = gsl_vector_get (x, 5); // normalized as muh0=Mh0

  // Matching scale
  long double mu = poi->Mt();
  
  MSinput MSm = MSinput::fromConsts( pow(mu,2), // Input scale
                                     muh,
                                     lam,
                                     yb, yt, g2, g1);
  
  // std::cout << "Mb= " << MSm.mb() << std::endl; 
  // std::cout << "MW= " << MSm.mW() << std::endl; 
  // std::cout << "MZ= " << MSm.mZ() << std::endl; 
  // std::cout << "MH= " << MSm.mH() << std::endl; 
  // std::cout << "Mt= " << MSm.mt() << std::endl;
  // std::cout << "V=  " << MSm.vev() << std::endl; 
  // std::cout << "1/al= " << 1./MSm.alpha() << std::endl; 


  long double a4pi = MSm.alpha()/4./Pi;
  long double as4pi = pow(gs/4./Pi,2);

  // b-quark mass
  // bb<MS> bQuark(MSm, pow(mu,2));
  // const double dMb = Mb0 - MSm.mb()*(1 +
  //                                    // EW
  //                                    a4pi*bQuark.x10() +
  //                                    a4pi*as4pi*bQuark.x11()+
  //                                    pow(a4pi,2)*bQuark.x20() +
  //                                    // QCD
  //                                    pow(as4pi,1)*bQuark.x01() +
  //                                    pow(as4pi,2)*bQuark.x02() +
  //                                    pow(as4pi,3)*bQuark.x03()
  //                                    );
  gsl_vector_set (f, 0, 0// dMb
                  );

  // W-boson mass
  WW<MS> WBoson(MSm, pow(mu,2));
  const double dMW = MW0 - MSm.mW()*(1 +
                                     // EW
                                     a4pi*WBoson.x10() +
                                     a4pi*as4pi*WBoson.x11()+
                                     pow(a4pi,2)*WBoson.x20()
                                     );
  gsl_vector_set (f, 1, dMW);

  // Z-boson mass
  ZZ<MS> ZBoson(MSm, pow(mu,2));
  const double dMZ = MZ0 - MSm.mZ()*(1 +
                                     // EW
                                     a4pi*ZBoson.x10() +
                                     a4pi*as4pi*ZBoson.x11()+
                                     pow(a4pi,2)*ZBoson.x20()
                                     );
  gsl_vector_set (f, 2, dMZ);

  // H-boson mass
  HH<MS> HBoson(MSm, pow(mu,2));
  const double dMH = MH0 - MSm.mH()*(1 +
                                     // EW
                                     a4pi*HBoson.x10() +
                                     a4pi*as4pi*HBoson.x11()+
                                     pow(a4pi,2)*HBoson.x20()
                                     );
  gsl_vector_set (f, 3, dMH);

  // t-quark mass
  tt<MS> tQuark(MSm, pow(mu,2));
  const double dMt = Mt0 - MSm.mt()*(1 +
                                     // EW
                                     a4pi*tQuark.x10() +
                                     a4pi*as4pi*tQuark.x11()+
                                     pow(a4pi,2)*tQuark.x20() +
                                     // QCD
                                     pow(as4pi,1)*tQuark.x01() +
                                     pow(as4pi,2)*tQuark.x02() +
                                     pow(as4pi,3)*tQuark.x03()
                                     );
  gsl_vector_set (f, 4, dMt);

  return GSL_SUCCESS;
}


int print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f % .3f % .3f % .3f % .3f "
          "f(x) = % .3e % .3e % .3e % .3e % .3e\n",
          iter,
          // 
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->x, 4),
          gsl_vector_get (s->x, 5),
          // 
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2),
          gsl_vector_get (s->f, 3),
          gsl_vector_get (s->f, 4),
          gsl_vector_get (s->f, 5));
}

SolvePMS::SolvePMS(OSinput oi)
{

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
  
  int status;
  size_t i, iter = 0;
  
  const size_t n = 6;

  gsl_multiroot_function f = {&Mdiff, n, &oi};


  double x_init[6] = {0.35761, 0.64822, // 1.1666, 
                      0.93558, 0, 0.12711, 132.03};
  gsl_vector *x = gsl_vector_alloc (n);

  for(int i = 0; i < n; i++)
    gsl_vector_set (x, i, x_init[i]);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 6);
  gsl_multiroot_fsolver_set (s, &f, x);

  print_state (iter, s);



  gsl_vector *fout = gsl_vector_alloc (n);

  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      
      print_state (iter, s);
      
      if (status)   /* check if solver is stuck */
        break;
      
      status =
        gsl_multiroot_test_residual (s->f, 1e-3);
    }
  while (status == GSL_CONTINUE && iter < 1000);
  
  printf ("status = %s\n", gsl_strerror (status));
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  
  
  // Mdiff(x, &oi, fout);

  // std::cout << "delta_MW_pole = " << gsl_vector_get (fout, 1) << std::endl;
  // std::cout << "delta_MZ_pole = " << gsl_vector_get (fout, 2) << std::endl;
  // std::cout << "delta_MH_pole = " << gsl_vector_get (fout, 3) << std::endl;
  // std::cout << "delta_Mt_pole = " << gsl_vector_get (fout, 4) << std::endl;
  
}

void SolvePMS::solve()
{

  double ms =90;


  // MSinput mi;
  
  // WW<MS> w(mi, ms);
  // ZZ<MS> z(mi, ms);
  // HH<MS> h(mi, ms);
  // tt<MS> t(mi, ms);
  
}



// RunUpto::RunUpto(OSinput oi, long double al_, long double as_, long double mu_): al0(al_), as0(as_), mu0(mu_)
// {
//   if(mu0 <  oi.Mb())
//     throw std::logic_error("ERROR: input at scale mu > Mb needed");

//   // We use matching at mu=Mt
//   // And need to run QCD and EW
//   // couplings up to Mt first

//   AlphaS       as;
//   AlphaQEDQCD  al;

//   // matching scale at top mass, when full SM aplicable
//   ms = oi.MMt();

//   WW<OS> w(oi, ms);
//   ZZ<OS> z(oi, ms);
//   HH<OS> h(oi, ms);
//   tt<OS> t(oi, ms);

//   long double aQCD = as(ms)/4./Pi;
//   long double aEW  = al.QED(ms)/4./Pi;

//   // 2-loop EW corrections
//   long double dWplus1 = 1 + aEW*w.y10() + aEW*aQCD*w.y11() + aEW*aEW*w.y20();
//   long double dZplus1 = 1 + aEW*z.y10() + aEW*aQCD*z.y11() + aEW*aEW*z.y20();
//   long double dHplus1 = 1 + aEW*h.y10() + aEW*aQCD*h.y11() + aEW*aEW*h.y20();
//   long double dtplus1 = 1 + aEW*t.y10() + aEW*aQCD*t.y11() + aEW*aEW*t.y20()
//     + aQCD*t.y01()+ aQCD*aQCD*t.y02()+ aQCD*aQCD*aQCD*t.y03();

//   long double Gf = pdg2014::Gf;

//   long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
//   long double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
    
//   a1 = 5./3.*(gg_ggp - gg)/16./Pi/Pi;
//   a2 = gg/16./Pi/Pi;
//   aS = aQCD;
//   ayt = pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
//   alam = Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;


//   std::cout << " At matching scale mu = " << ms << std::endl;
//   std::cout << " g1 = " << sqrt(3./5.*a1)*4*Pi << std::endl;
//   std::cout << " g2 = " << sqrt(a2)*4*Pi << std::endl;
//   std::cout << " g3 = " << sqrt(aS)*4*Pi << std::endl;
//   std::cout << " yt = " << sqrt(ayt)*4*Pi << std::endl;
//   std::cout << " lam = " << alam*16*Pi*Pi << std::endl;

//   av = new CouplingsSM<3,3,3,3,-1,-1,3>(a1,a2,aS,ayt,0,0,alam,pow(ms,2),3);
//   // std::cout << " LAMMMMMMMMMMMM " << av->operator()(pow(10000,2))[6] << " ::: " << lambda(10000) << std::endl;
// }

// state_type RunUpto::operator()(long double mu)
// {  
//   return av->operator()(pow(mu,2));
// }


// long double RunUpto::lambda(long double mu)
// {
//   return av->operator()(pow(mu,2))[6];
// }
