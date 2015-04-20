#include <iostream>
#include <cmath>
#include "mr.hpp"
#include "tools.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Default log level is logERROR
      // loglevel = logDEBUG;
            
      // Input: Pole masses and Fermi constant in OS scheme
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);

      // Running QCD coupling for as(Mt) from as(MZ)
      AlphaS as(oi);

      // Set of all running parameters at scale Mt
      P2MS pMSmt(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::all);

      // Initial values input by hand
      Couplings<3,3,3,
                3,3,-1,
                3,3,0> av(
                          5./3.*pow(0.35830/4./Pi,2), // GUT normalization
                          pow(0.64779/4./Pi,2),
                          pow(1.1666/4./Pi,2),
                          pow(0.93690/4./Pi,2),
                          pow(0.0/4./Pi,2),
                          pow(0.0/4./Pi,2),
                          0.12604*pow(4.*Pi,-2),
                          131.55,
                          0*246,
                          pow(173.34,2)
                          );
      
      Couplings<3,3,3,
                3,0,-1,
                3,3,0> avP2MS(pMSmt);
      
      
      lout(logINFO) << "Prepared for plotting!";
      std::ofstream fout("cEvol.dat");
      
      // Header
      fout << "# scale a1 a2 a3 at ab alam mu0/1000\n";
      for (double muPow = 2.; muPow <= 20.; muPow+=0.5)
        {

          lout(logINFO) << "Scale is 10^" << muPow;

          state_type v = avP2MS(pow(10,2*muPow));
          
          fout << muPow << " "
               << sqrt(v[0])*4.*Pi          << " "
               << sqrt(v[1])*4.*Pi          << " "
               << sqrt(v[2])*4.*Pi          << " "
               << sqrt(v[3])*4.*Pi          << " "
               << sqrt(v[4])*4.*Pi          << " "
            // << v[5]          << " "
               << v[6]*16.*Pi*Pi            << " "
               << v[7]/1000. << std::endl;

        }
      
      
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

