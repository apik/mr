#include <iostream>
#include <cmath>
#include <fstream>
#include "mr.hpp"
#include "tools.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Default log level is logERROR
      loglevel = logDEBUG;
#define exType 2


#if exType==1
      // Input: Pole masses and Fermi constant in OS scheme
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);

      // Running QCD coupling for as(Mt) from as(MZ)
      AlphaS as(oi);

      // Set of all running parameters at scale Mt
      P2MS pMSmt(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::all);

      std::cout << "Critical Top mass at lambda(Mpl)=0 is " << critMt_scaleNotFixed(oi, pdg2014::Mpl) << std::endl;
        
#elif exType==2

      std::ofstream fout("cMHMtCrit.dat");

      // Start point
      std::pair<long double,long double> metaStable_Stable(pdg2014::Mt, pdg2014::Mpl);
      std::pair<long double,long double> metaStable_Instable(pdg2014::Mt, pdg2014::Mpl);
      for(int MHPoint = 120; MHPoint <= 132; MHPoint += 2)
        {
          OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, double(MHPoint), metaStable_Stable.first);

          // alpha_S error bars
          metaStable_Stable = critMt_scaleNotFixed<Mt_Stability>(oi, metaStable_Stable.second, pdg2014::asMZ);
          std::pair<long double,long double> metaStable_StableL = critMt_scaleNotFixed<Mt_Stability>(oi, metaStable_Stable.second, pdg2014::asMZ - 0.0006);
          std::pair<long double,long double> metaStable_StableU = critMt_scaleNotFixed<Mt_Stability>(oi, metaStable_Stable.second, pdg2014::asMZ + 0.0006);
          
          // metaStable_Instable = 0;// critMt_scaleNotFixed<Mt_Instability>(oi, pdg2014::Mpl);

          std::cout << "Critical Top mass    at Mh = " << oi.MH() << " is " << metaStable_Stable.first   << std::endl; 
          // std::cout << "Instability Top mass at Mh = " << oi.MH() << " is " << metaStable_Instable << std::endl; 
          fout << oi.MH() << " " << metaStable_StableL.first << " " << metaStable_Stable.first<< " " << metaStable_StableU.first << " " << metaStable_Stable.second << " " << metaStable_Instable.first << std::endl;
          
        }
#endif

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

