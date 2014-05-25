#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
 
      // Mihaila, Salomon, Steinhauser
      CouplingsSM<3,3,3,3,3,3,3> a333(
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
      


      Plot2 plotMratio("mratio", "mH/mW, mH/mt", "log_10(\\mu/GeV)", "mH/mW, mH/mt", "mH/mW", "mH/mt");
      // Input from Degrassi,Strumia,...
      SMinput DS2l(0, 80.384, 91.1876, 125.7, 173.10);
      for (int lg10mu = 2; lg10mu < 21; lg10mu++)
        {

          std::cout << " Scale \\mu = 10^" << lg10mu << " GeV:" << std::endl;

          WW dW  = WW(DS2l, pow(10.,2.*lg10mu));          
          HH dH  = HH(DS2l, pow(10.,2.*lg10mu));          
          tt dt  = tt(DS2l, pow(10.,2.*lg10mu));          

          // Couplings:
          state_type av = a333(pow(10.,2.*lg10mu));
          long double aEM = av[0]*av[1]/(av[0]+av[1]);
          long double aS  = av[2];
          std::cout << " a_EM = " << aEM << std::endl;
          std::cout << " a_yt = " << aEM << std::endl;
          std::cout << " a_lam = " << av[6] << std::endl;
          
          long double mmW = DS2l.MMW() * ( 
                                          1 +
                                          aEM *       dW.m10().real() +  
                                          aEM * aS  * dW.m11().real() +
                                          aEM * aEM * dW.m20().real() );

          std::cout << " mmW = " << mmW << std::endl;

          long double mmH = DS2l.MMH() * ( 
                                          1 +
                                          aEM *       dH.m10().real() +  
                                          aEM * aS  * dH.m11().real() +
                                          aEM * aEM * dH.m20().real() );

          std::cout << " mmH = " << mmH << std::endl;
 
          long double mt = DS2l.Mt() * ( 
                                        1 +
                                        aEM *       dt.m10().real() +  
                                        aEM * aS  * dt.m11().real() +
                                        aEM * aEM * dt.m20().real() );

          std::cout << " mt = " << mt << std::endl;
         
          // plotMratio.add(lg10mu,sqrt(6*mmH/mmW),sqrt(6*mmH)/mt);

          plotMratio.add(lg10mu, 
                         sgn(av[6]) * sqrt(4. * fabs(av[6])/av[1]),  // mH/mW
                         sgn(av[6]) * sqrt(8. * fabs(av[6])/av[3])); // mH/mt
        }


    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

