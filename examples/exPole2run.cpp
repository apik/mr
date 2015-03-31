#include <iostream>
#include <cmath>
#include "mr.hpp"
#include "tools.hpp"
#include "gnuplot.hpp"

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
                3,0,0> avP2MS(pMSmt);



      Couplings<1,1,3,
                3,3,-1,
                3,0,0> avP2MS1(pMSmt);
      
      Couplings<2,2,3,
                3,3,-1,
                3,0,0> avP2MS2(pMSmt);
      
      Couplings<3,3,3,
                3,3,-1,
                3,0,0> avP2MS3(pMSmt);
      
      state_type avMplP2MS4 = avP2MS(mmPlanck/4.);
      state_type avMplP2MS3 = avP2MS(mmPlanck/3.);
      state_type avMplP2MS2 = avP2MS(mmPlanck/2.);
      state_type avMplP2MS = avP2MS(mmPlanck);

      std::cout << "muP2MS = " << avMplP2MS[7] << std::endl;


      critMH(oi, mmPlanck);
        
      // Plot3 plotGL("runSM_gauge", "3-loop gauge couplings running", "log_10(\\mu/GeV)", "a1,a2,a3", "a1", "a2", "a3");
      
      // for (double  lg10mu = 17; lg10mu > 2; lg10mu -= 0.1)
      //   {
          
      //     std::cout << "Scale: " << powf(10.,2*lg10mu) << "         And Planck= " << mmPlanck << std::endl;
          
      //     // std::cout << "               be1 = " << be1(aSM0,aSM,1000.) << std::endl;
          
      //     state_type avFull = avP2MS(powf(10.,2*lg10mu));
          
      //     // std::cout << "a1 = " << 4.*Pi*sqrt(avFull[0]) << std::endl;
      //     // std::cout << "a2 = " << 4.*Pi*sqrt(avFull[1]) << std::endl;
      //     // std::cout << "a3 = " << 4.*Pi*sqrt(avFull[2]) << std::endl;
      //     // std::cout << "at = " << 4.*Pi*sqrt(avFull[3]) << std::endl;
      //     // std::cout << "lm = " << 16.*Pi*Pi*avFull[6] << std::endl;
      //     // std::cout << "mu = " << avFull[7] << std::endl;
      //     // std::cout << "v = " << avFull[8] << std::endl << std::endl;

      //     std::cout << "a1 = " << 4.*Pi*sqrt(avMplP2MS[0]) << std::endl;
      //     std::cout << "a2 = " << 4.*Pi*sqrt(avMplP2MS[1]) << std::endl;
      //     std::cout << "a3 = " << 4.*Pi*sqrt(avMplP2MS[2]) << std::endl;
      //     std::cout << "at = " << 4.*Pi*sqrt(avMplP2MS[3]) << std::endl;
      //     std::cout << "lm = " << 16.*Pi*Pi*avMplP2MS[6] << std::endl;
      //     std::cout << "mu = " << avMplP2MS[7] << std::endl;
      
      //     std::cout << "v = " << avMplP2MS[8] << std::endl << std::endl;

          
      //     plotGL.add(lg10mu,4.*Pi*avFull[0], 4.*Pi*avFull[1], 4.*Pi*avFull[2]);
      //     // std::cout << "Point: Log[mu] = " << lg10mu << " lam = " << 16.*Pi*Pi*avFull[6] << std::endl;
      //   }


      // std::string legends[6] = {"a1 1-loop", "a1 2-loop", "a1 3-loop", "a2 1-loop", "a2 2-loop", "a2 3-loop"};

      // Plot<6> plotWeyl("runSM_gaugeZoom", "\\lambda running in 3-3-3, 3-3-2 and 3-2-1 schemes", "log_10(\\mu/GeV)", "a1,a2",legends );

      
      // for (double  lg10mu = 13.; lg10mu < 13.1; lg10mu += 0.01)
      //   {
          
      //     std::cout << "Scale: " << powf(10.,2*lg10mu) << "         And Planck= " << mmPlanck << std::endl;
          
      //     state_type av1 = avP2MS1(powf(10.,2*lg10mu));
      //     state_type av2 = avP2MS2(powf(10.,2*lg10mu));
      //     state_type av3 = avP2MS3(powf(10.,2*lg10mu));
          

      //     long double pltv[6] = {4.*Pi*av1[0],
      //                            4.*Pi*av2[0],
      //                            4.*Pi*av3[0],
      //                            4.*Pi*av1[1],
      //                            4.*Pi*av2[1],
      //                            4.*Pi*av3[1]};
          
      //     plotWeyl.add(lg10mu,pltv);
      //     // std::cout << "Point: Log[mu] = " << lg10mu << " lam = " << 16.*Pi*Pi*avFull[6] << std::endl;
      //   }


      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

