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
      typedef CouplingsSM<3,3,3,3,3,3,3> QCD;
      typedef CouplingsMu QCD1;



      //// 
      ///   3-loop Gauge couplings running
      // 

      // Mihaila, Salomon, Steinhauser
      QCD1 qcd1_(
                  0.,
                  0.,
                  0.108057/4./Pi,
                  0.,
                  0.,
                  0.,
                  0.013/pow(4.*Pi,2),             // Lambda
		  0, //mu2
                  pow(172.1,2),
                  3 // NG
                  );
      
      QCD qcd_(
                  0.,
                  0.,
                  0.108057/4./Pi,
                  0.,
                  0.,
                  0.,
                  0.013/pow(4.*Pi,2),             // Lambda
                  pow(172.1,2),
                  3 // NG
	      );
      // 
      
      state_type qcd_mz = qcd_(pow(91.19,2));
      state_type qcd1_mz = qcd1_(pow(91.19,2));

      std::cout << "as6(mz)" << (4.*Pi)*qcd_mz[2] << std::endl;
      std::cout << "as6(mz)" << (4.*Pi)*qcd1_mz[2] << std::endl;
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

