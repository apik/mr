#include <iostream>
#include <cmath>
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
      P2MS pMSmt(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::all);

      // Set of all running parameters at scale MZ
      P2MS pMSmZ(oi,pdg2014::Gf, as(oi.MZ()), oi.MZ(), order::all);

      
      // Input at mu=MZ for running as in [hep-ph]1208.3357
      // std::cout << "alpha(1) = " << parsMS.a1()*4*Pi << std::endl
      //           << "alpha(2) = " << parsMS.a2()*4*Pi << std::endl
      //           << "alpha(3) = " << parsMS.as()*4*Pi << std::endl
      //           << "alpha(t) = " << parsMS.at()*4*Pi << std::endl
      //           << "alpha(b) = " << parsMS.ab()*4*Pi << std::endl
      //           << "4*Pi*lam = " << parsMS.alam()*pow(4*Pi,2) << std::endl;





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
      

      long double mmPlanck = pow(1.2209,2) * pow(10.,2*19); 
      state_type avMpl = av(mmPlanck);

      // std::cout << avMpl[8]

      std::cout << "a1 = " << 4.*Pi*sqrt(avMpl[0]) << std::endl;
      std::cout << "a2 = " << 4.*Pi*sqrt(avMpl[1]) << std::endl;
      std::cout << "a3 = " << 4.*Pi*sqrt(avMpl[2]) << std::endl;
      std::cout << "at = " << 4.*Pi*sqrt(avMpl[3]) << std::endl;
      std::cout << "lm = " << 16.*Pi*Pi*avMpl[6] << std::endl;
      std::cout << "mu = " << avMpl[7] << std::endl;
      
      std::cout << "v = " << avMpl[8] << std::endl;
      std::cout << "vc= " << avMpl[7]/sqrt(2.*16.*Pi*Pi*avMpl[6]) << std::endl;



      Couplings<3,3,3,
                3,3,-1,
                3,3,0> avP2MS(pMSmt);

      state_type avMplP2MS = avP2MS(mmPlanck);

      std::cout << "muP2MS = " << avMplP2MS[7] << std::endl;
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

