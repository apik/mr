#include <iostream>
#include <cmath>
#include "mr.hpp"

using namespace mr;

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
      
      // Initial values for running, input from pole masses
      Couplings<3,3,3,
        3,0,-1,
        3,3,0> avP2MS(pMSmt);
      
      std::cout << std::setprecision(3);
      
      for (size_t muPow = 3.; muPow <= 20.; muPow++)
        {

          SMCouplings av = avP2MS(pow(10,2*muPow));
          
          std::cout << " log10(mu) = " << muPow 
                    << " a1 = " << av[couplings::g1]
                    << " a2 = " << av[couplings::g2]
                    << " a3 = " << av[couplings::gs]
                    << " at = " << av[couplings::yt]
                    << " ab = " << av[couplings::yb]
                    << " atau = " << av[couplings::ytau]
                    << " alam = " << av[couplings::lam]
                    << " mu0 = " << av[couplings::mu0] << std::endl;
        }      
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

