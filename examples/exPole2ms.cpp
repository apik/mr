// Example of Pole masses and Gf conversion to
// set of running couplings, running Higgs mass
// term and running vev in MS scheme

#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {

      // Input: Pole masses and Fermi constant in OS scheme
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);

      // Running QCD coupling for as(Mt) from as(MZ)
      AlphaS as(oi);

      // Set of all running parameters at scale Mt
      P2MS parsMS(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::all);

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

