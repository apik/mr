#include <iostream>
#include <cmath>
#include "mr.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {

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


      ////
      ///  QCD + Yukawa top + lambda
      //

      SM_Ht aHt(
                   0,
                   0,
                   pow(1.1646028932952313/4./Pi,2),
                   pow(0.9344417846242093/4/Pi,2),
                   0,
                   0,
                   0.13649339025452764*pow(4.*Pi,-2),
                   pow(172.9,2),
                   3);


    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

