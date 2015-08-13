#include <iostream>
#include <cmath>
#include "mr.hpp"

using namespace mr;

int main (int argc, char *argv[])
{
  try
    {
      // Default log level is logERROR
      loglevel = logDEBUG;

      // Input: Pole masses and Fermi constant in OS scheme
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);

      // Running QCD coupling for as(Mt) from as(MZ)
      AlphaS as(oi);

      // Set of all running parameters at scale Mt
      P2MS pMSmt(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::all);

      std::cout << "Critical Higgs mass at lambda(Mpl)=0 is " << critMH_scaleNotFixed(oi, pdg2014::Mpl) << std::endl;

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

