#include <iostream>
#include <cmath>
#include "mr.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {

      // Jegerlehner test:

      // MZ
      // CouplingsVevMu av(
      //                   5./3.*pow(0.3497/4./Pi,2), // GUT normalization
      //                   pow(0.6530/4./Pi,2),
      //                   pow(1.2200/4./Pi,2),
      //                   pow(0.9347/4./Pi,2),
      //                   pow(0.0238/4./Pi,2),
      //                   pow(0.0104/4./Pi,2),
      //                   0.8070/6.*pow(4.*Pi,-2),
      //                   sqrt(3./0.8070)*89.096*sqrt(2.),
      //                   89.096*sqrt(2.),
      //                   pow(91.1876,2),
      //                   3);

      // Mt Degrassi
      CouplingsVevMu av(
                        5./3.*pow(0.3587/4./Pi,2), // GUT normalization
                        pow(0.6483/4./Pi,2),
                        pow(1.1644/4./Pi,2),
                        pow(0.9399/4./Pi,2),
                        0*pow(0.0238/4./Pi,2),
                        0*pow(0.0104/4./Pi,2),
                        0.7626/6.*pow(4.*Pi,-2),
                        0*sqrt(3./0.7626)*97.278*sqrt(2.),
                        97.278*sqrt(2.),
                        pow(173.5,2),
                        3);

      
      // Mt
      // CouplingsVevMu av(
      //                   5./3.*pow(0.3509/4./Pi,2), // GUT normalization
      //                   pow(0.6496/4./Pi,2),
      //                   pow(1.1644/4./Pi,2),
      //                   pow(0.9002/4./Pi,2),
      //                   pow(0.0227/4./Pi,2),
      //                   pow(0.0104/4./Pi,2),
      //                   0.7373/6.*pow(4.*Pi,-2),
      //                   sqrt(3./0.7373)*89.889*sqrt(2.),
      //                   89.889*sqrt(2.),
      //                   pow(173.5,2),
      //                   3);


      long double mmt = pow(173.5,2); 
      state_type  amt = av(mmt);

      std:: cout << "g1 at mmt = " << 4.*Pi*sqrt(3/5.*amt[0]) << std::endl;
      std:: cout << "g2        = " << 4.*Pi*sqrt(amt[1]) << std::endl;
      std:: cout << "g3        = " << 4.*Pi*sqrt(amt[2]) << std::endl;
      std:: cout << "yt        = " << 4.*Pi*sqrt(amt[3]) << std::endl;
      std:: cout << "yb        = " << 4.*Pi*sqrt(amt[4]) << std::endl;
      std:: cout << "ytau      = " << 4.*Pi*sqrt(amt[5]) << std::endl;
      std:: cout << "lamda     = " << 6*pow(4.*Pi,2)*amt[6] << std::endl;
      std:: cout << "v         = " << amt[7] << std::endl; 
      std:: cout << "mu        = " << amt[8]/sqrt(2) << std::endl; 
      std:: cout << "------------" << std::endl;
      std:: cout << "calc v    = " << sqrt(3./(6*pow(4.*Pi,2)*amt[6]))*amt[8] << std::endl;
      std:: cout << std::endl;

      
      long double mu0 = pow(1.1,2) * pow(10.,2*16); 
      state_type  amu0 = av(mu0);

      std:: cout << "g1 at mu0 = " << 4.*Pi*sqrt(3/5.*amu0[0]) << std::endl;
      std:: cout << "g2        = " << 4.*Pi*sqrt(amu0[1]) << std::endl;
      std:: cout << "g3        = " << 4.*Pi*sqrt(amu0[2]) << std::endl;
      std:: cout << "yt        = " << 4.*Pi*sqrt(amu0[3]) << std::endl;
      std:: cout << "yb        = " << 4.*Pi*sqrt(amu0[4]) << std::endl;
      std:: cout << "ytau      = " << 4.*Pi*sqrt(amu0[5]) << std::endl;
      std:: cout << "lamda     = " << 6*pow(4.*Pi,2)*amu0[6] << std::endl;
      std:: cout << "v         = " << amu0[7] << std::endl; 
      std:: cout << "mu        = " << amu0[8]/sqrt(2) << std::endl; 
      std:: cout << "------------" << std::endl;
      std:: cout << "calc v    = " << sqrt(3./(6*pow(4.*Pi,2)*amu0[6]))*amu0[8] << std::endl;
      std:: cout << std::endl;
      
      
      long double mmPlanck = pow(1.2209,2) * pow(10.,2*19); 
      state_type  aMpl = av(mmPlanck);

      std:: cout << "g1 at Mpl = " << 4.*Pi*sqrt(3/5.*aMpl[0]) << std::endl;
      std:: cout << "g2        = " << 4.*Pi*sqrt(aMpl[1]) << std::endl;
      std:: cout << "g3        = " << 4.*Pi*sqrt(aMpl[2]) << std::endl;
      std:: cout << "yt        = " << 4.*Pi*sqrt(aMpl[3]) << std::endl;
      std:: cout << "yb        = " << 4.*Pi*sqrt(aMpl[4]) << std::endl;
      std:: cout << "ytau      = " << 4.*Pi*sqrt(aMpl[5]) << std::endl;
      std:: cout << "lamda     = " << 6*pow(4.*Pi,2)*aMpl[6] << std::endl;
      std:: cout << "v         = " << aMpl[7] << std::endl; 
      std:: cout << "mu        = " << aMpl[8]/sqrt(2) << std::endl;
      std:: cout << "------------" << std::endl;
      std:: cout << "calc v    = " << sqrt(3./(6*pow(4.*Pi,2)*aMpl[6]))*aMpl[8] << std::endl;
      std:: cout << std::endl;

      return 0;
      
      // Standard Model without bottom and tau Yukawa Couplings, 
      // 5 active couplings used in [arXiv:1205.2893]
      // available in F.Bezrukov code from http://www.inr.ac.ru/~fedor/SM/
      typedef CouplingsSM<3,3,3,3,-1,-1,3> SM_5;


      // Mihaila, Salomon, Steinhauser (MSS), [arXiv:1208.3357]
      // Only gauge 3-loop beta-functions available
      typedef CouplingsSM<3,3,3,2,2,2,2> SM_MSS;


      // Running motivated by Weyl consistancy conditions
      // Instead of 3-3-3 running use 3-2-1, [arXiv:1306.3234]
      typedef CouplingsSM<3,3,3,2,-1,-1,1> SM_Weyl;


      // Leading 3-loop corrections in QCD with Higgs self-coupling
      // and Yukawa top, from  [arXiv:1205.2892]
      typedef CouplingsSM<-1,-1,3,3,-1,-1,3> SM_Ht;

      

      // Input from Chetyrkin, Zoller
      SM_5 aFull(
                 5./3.*pow(0.3587/4./Pi,2),
                 pow(0.6484/4./Pi,2),
                 pow(1.1644/4./Pi,2),
                 pow(0.94/4./Pi,2),
                 pow(0.02/4./Pi,2),
                 pow(0.01/4./Pi,2),
                 0.13*pow(4.*Pi,-2),
                 pow(173.5,2),
                 3);
      
      CouplingsSM<2,2,3,3,-1,-1,3> aHt(
                 5./3.*pow(0.3587/4./Pi,2),
                 pow(0.6484/4./Pi,2),
                 pow(1.1644/4./Pi,2),
                 pow(0.94/4./Pi,2),
                 pow(0.02/4./Pi,2),
                 pow(0.01/4./Pi,2),
                 0.13*pow(4.*Pi,-2),
                 pow(173.5,2),
                 3);

      CouplingsSM<2,2,2,2,2,2,2> a2(
                 5./3.*pow(0.3587/4./Pi,2),
                 pow(0.6484/4./Pi,2),
                 pow(1.1644/4./Pi,2),
                 pow(0.94/4./Pi,2),
                 pow(0.02/4./Pi,2),
                 pow(0.01/4./Pi,2),
                 0.13*pow(4.*Pi,-2),
                 pow(173.5,2),
                 3);

          // Check beta

//  b1 = -3.292345e-01

//    Final model parameters:
// Q      = 15847915711.7119034147
// m2     =  0.0000000000
// lambda = -0.0000267959
// yt     =  0.5360697333
// g3     =  0.6470656003
// g      =  0.5680742232
// gp     =  0.4026628421
// v      =  0.0000000000

          state_type aSM0,betaSM(7);
          // aSM0.push_back(0.4026628421);
          // aSM0.push_back(0.5680742232);
          // aSM0.push_back(0.6470656003);
          // aSM0.push_back(0.5360697333);
          // aSM0.push_back(0);
          // aSM0.push_back(0);
          // aSM0.push_back(-0.0000267959);


          aSM0.push_back(5./3.*pow(0.3587/4./Pi,2));
          aSM0.push_back(pow(0.6484/4./Pi,2));
          aSM0.push_back(pow(1.1644/4./Pi,2));
          aSM0.push_back(pow(0.94/4./Pi,2));
          aSM0.push_back(0);
          aSM0.push_back(0);
          aSM0.push_back(0.13*pow(4.*Pi,-2));



          BetaSMFull(1, 1, 1, 1, -1, -1, 1, 3)(aSM0,betaSM,1000.);
          std::cout << "               be1 = " << 3./5.*betaSM[0] << " be2 = " << betaSM[1] << " be3 = " << betaSM[2] << std::endl;
          std::cout << "               be7 = " << 16.*Pi*Pi*betaSM[6]*2. << " be4 = " << betaSM[3]  << std::endl;

          // BetaSMFull(-1, -1, -1, -1, -1, -1, 2, 3)(aSM0,betaSM,1000.);
          // std::cout << "               be2 = " << betaSM[6] << std::endl;
          // BetaSMFull(-2, -2, -2, -2, -2, -2, 3, 3)(aSM0,betaSM,1000.);
          // std::cout << "               be3 = " << betaSM[6] << std::endl;

          // return 0;


      Plot3 plotGL("runSM_MSS_gauge", "3-loop gauge couplings running", "log_10(\\mu/GeV)", "lam", "lambda full", "lambda H,t", "lambda 2-loop");

      for (double  lg10mu = 10.32; lg10mu < 10.33; lg10mu += 0.00001)
        {


          
          // std::cout << "               be1 = " << be1(aSM0,aSM,1000.) << std::endl;
          
          state_type avFull = aFull(powf(10.,2*lg10mu));

          state_type avHt = aHt(powf(10.,2*lg10mu));

          state_type av2 = a2(powf(10.,2*lg10mu));
          
          plotGL.add(lg10mu,16.*Pi*Pi*avFull[6], 16.*Pi*Pi*avHt[6], 16.*Pi*Pi*av2[6]);
          std::cout << "Point: Log[mu] = " << lg10mu << " lam = " << 16.*Pi*Pi*avFull[6] << std::endl;
        }

      return 0;
      

      //// 
      ///   3-loop Gauge couplings running
      // 

      // Mihaila, Salomon, Steinhauser
      SM_MSS aMSS(
                  0.0169225/4./Pi,
                  0.033735/4./Pi,
                  0.1173/4./Pi,
                  0.07514/4./Pi,
                  0.00002064/4./Pi,
                  8.077E-6/4./Pi,
                  0.13/pow(4.*Pi,2),             // Lambda
                  pow(91.1876,2),
                  3 // NG
                  );
      
      Plot3 plotMSS("runSM_MSS_gauge", "3-loop gauge couplings running", "log_10(\\mu/GeV)", "a1,a2,a3", "a1 ", "a2", "a3");

      for (int lg10mu = 1; lg10mu < 18; lg10mu++)
        {
          state_type av = aMSS(pow(10.,2.*lg10mu));
          
          plotMSS.add(lg10mu,4.*Pi*av[0],4.*Pi*av[1],4.*Pi*av[2]);
        }

      /////
      ////  Compare 3-3-3 and 3-2-1 running
      ///   
      // 
      
      std::cout << "Beta functions for 3-3-3 running" << std::endl;
      // Weyl
      SM_5 a333(
                5./3.*pow(0.358729/4./Pi,2),
                pow(0.648382/4/Pi,2),
                pow(1.1646/4/Pi,2),
                pow(0.934442/4/Pi,2),
                0,
                0,
                0.136493/pow(4.*Pi,2),             // Lambda
                pow(172.9,2),
                3 // NG
                );
      std::cout << "Beta functions for 3-2-1 running" << std::endl;
      SM_Weyl a321(
                   5./3.*pow(0.358729/4./Pi,2),
                   pow(0.648382/4/Pi,2),
                   pow(1.1646/4/Pi,2),
                   pow(0.934442/4/Pi,2),
                   0,
                   0,
                   0.136493/pow(4.*Pi,2),             // Lambda
                   pow(172.9,2),
                   3 // NG
                   );
      
      Plot2 plotWeyl("runSM_Weyl", "\\lambda running in 3-3-3 and 3-2-1 schemes", "log_10(\\mu/GeV)", "\\lambda", "3-3-3 ", "3-2-1");

      for (int lg10mu = 8; lg10mu < 19; lg10mu++)
        {
          state_type av333 = a333(pow(10.,2.*lg10mu));
          state_type av321 = a321(pow(10.,2.*lg10mu));
          
          plotWeyl.add(lg10mu,pow(4.*Pi,2)*av333[6],pow(4.*Pi,2)*av321[6]);
        }


      ////
      /// Jegerlehner analysis [arXiv:1304.7813]
      // 
      SM_5 aSM_FJ(5/3*pow(0.35/4/Pi,2),
                  pow(0.653/4/Pi,2),
                  pow(1.22/4/Pi,2),
                  pow(0.935/4/Pi,2),
                  0,
                  0,
                  0.132667/pow(4.*Pi,2),
                  pow(91.1876,2), // Initial Scale MZ
                  3
                  );      


      // ////
      // ///  QCD + Yukawa top + lambda
      // //

      // SM_Ht aHt(
      //              0,
      //              0,
      //              pow(1.1646028932952313/4./Pi,2),
      //              pow(0.9344417846242093/4/Pi,2),
      //              0,
      //              0,
      //              0.13649339025452764*pow(4.*Pi,-2),
      //              pow(172.9,2),
      //              3);


    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

