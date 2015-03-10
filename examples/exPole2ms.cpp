#include <iostream>
#include <cmath>
#include "mr.hpp"


int main (int argc, char *argv[])
{
  try
    {

      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);
      AlphaS as(oi);
      
      P2MS spms(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::allQCD);




      std::cout << "Orders all: " << std::hex << (order::x01|order::x10|order::x02|order::x11|order::x20|order::x03) <<std::endl;
      
      
      std::cout << "Orders QCD all: " << std::hex << (order::x01|order::x02|order::x03) <<std::endl;
      
      std::cout << "Orders all EW: " << std::hex << (order::x10|order::x11|order::x20) <<std::endl;
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

