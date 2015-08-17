#include <iostream>
#include <fstream>
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
      
      std::cout << "Critical Higgs mass from condition lam=0 at Planck scale      " 
                << critMH(oi, pdg2014::Mpl) << std::endl;

      std::cout << "Critical Higgs mass from condition beta_lam=0 at Planck scale " 
                << critMH(oi, pdg2014::Mpl) << std::endl;


    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

