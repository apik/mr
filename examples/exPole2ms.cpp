#include <iostream>
#include <cmath>
#include "mr.hpp"

#include "p2ms.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {

      OSinput oi(0, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);
      AlphaS as(oi);
      
      P2MS spms(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt());
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

