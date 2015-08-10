#include <iostream>
#include <cmath>
// #include <boost/numeric/odeint.hpp>

#include "mr.hpp"
#include "CRunDec.h"

int main (int argc, char *argv[])
{
  try
    {

      OSinput oi(pdg2014::mb,pdg2014::MW,pdg2014::MZ,pdg2014::MH,pdg2014::Mt);
      
      AlphaS aoi(oi);

      

      
      long double alphaMt  = 0.00779305;
      long double alphaS = 0.1184;

      // Default loops=4, nf=5
      AlphaS as4l5(4,5);
      AlphaS as3l5(3,5);
      AlphaS as2l5(2,5);
      AlphaS as1l5(1,5);

      // CRunDec for comparison nf=5
      CRunDec* pObjnf5 = new CRunDec(5);      

      std::cout.precision(20);
      std::cout << std::fixed;

      std::cout << "loops = 4, nf = 5, \\alpha_s(MMt) = " << as4l5(pow(pdg2014::Mt,2)) << std::endl;

      long double mu20 = pow(pdg2014::MZ,2);
      long double step = (pow(pdg2014::Mt,2) - pow(pdg2014::MZ,2))/10;
      
      // for(int i = 1; i < 10; i++)
      //   {
      //     long double mu2 = mu20 + i*step;
      
      //     std::cout << "Difference for \\mu=" << sqrt(mu2) << std::endl 
      //               << "\tin 4-loop: " << as4l5(mu2) - pObjnf5 -> AlphasExact(alphaS,pdg2014::MZ,sqrt(mu2),4) << " "
      //               << 100*fabs((as4l5(mu2) - pObjnf5 -> AlphasExact(alphaS,pdg2014::MZ,sqrt(mu2),4))/as4l5(mu2)) << " % " << std::endl

      //               << "\tin 3-loop: " << as3l5(mu2) - pObjnf5 -> AlphasExact(alphaS,pdg2014::MZ,sqrt(mu2),3) << " "
      //               << 100*fabs((as3l5(mu2) - pObjnf5 -> AlphasExact(alphaS,pdg2014::MZ,sqrt(mu2),3))/as3l5(mu2)) << " % "<< std::endl

      //               << "\tin 2-loop: " << as2l5(mu2) - pObjnf5 -> AlphasExact(alphaS,pdg2014::MZ,sqrt(mu2),2) << " "
      //               << 100*fabs((as2l5(mu2) - pObjnf5 -> AlphasExact(alphaS,pdg2014::MZ,sqrt(mu2),2))/as2l5(mu2)) << " % "<< std::endl

      //               << "\tin 1-loop: " << as1l5(mu2) - pObjnf5 -> AlphasExact(alphaS,pdg2014::MZ,sqrt(mu2),1) << " "
      //               << 100*fabs((as1l5(mu2) - pObjnf5 -> AlphasExact(alphaS,pdg2014::MZ,sqrt(mu2),1))/as1l5(mu2)) << " % "<< std::endl << std::endl;


      //     std::cout << " My : " << as4l5(mu2) << std::endl; 
      //     std::cout << " Aoi: " << aoi.run(0.1184, pdg2014::MZ, sqrt(mu2), 5) << std::endl; 
      //   }


      
      // as4l5.upto(3);

      // std::cout << as4l5.upto(3) << std::endl;

      CRunDec crundec(5);
      crundec.nfMmu[0].nf = 6;
      crundec.nfMmu[0].Mth = pdg2014::Mt;
      crundec.nfMmu[0].muth = pdg2014::Mt;
      
      crundec.AlL2AlH(pdg2014::asMZ, pdg2014::MZ,crundec.nfMmu,pdg2014::Mt,4);
      
      std::cout << "My code: \\mu=M_top, g_3 = " << as4l5(pow(pdg2014::Mt,2)) << std::endl;
      std::cout << "CRunDec: \\mu=M_top, g_3 = " << crundec.AlL2AlH(pdg2014::asMZ, pdg2014::MZ,crundec.nfMmu,pdg2014::Mt,4) << std::endl;

      std::cout << " Aoi(50): " << run(pdg2014::asMZ, pdg2014::MZ, 50., 5, 4 ) << std::endl; 

      std::cout << " Dec: " << 4.*Pi*as5nf2as6nf(pow(pdg2014::Mt,2), pow(100*pdg2014::Mt,2), pdg2014::asMZ/4./Pi, /* nl= */5, 3) << std::endl;


      long double aMt5,aMt6,ak6;


      std::cout <<  "1-loop decoupling " <<std::endl;
      aMt5 = run(pdg2014::asMZ, pdg2014::MZ, pdg2014::Mt, 5 , 4);
      
      aMt6 = as5nf2as6nf(pdg2014::Mt, pdg2014::Mt, aMt5, /* nl= */5, 1);
      
      ak6 = run(aMt6, pdg2014::Mt, 100000, 6 );

      std::cout << " at5 = " << aMt5 << "  at6 = " << aMt6  << "  ak6 = " << ak6 << std::endl;

      std::cout <<  "2-loop decoupling " <<std::endl;
      aMt5 = run(pdg2014::asMZ, pdg2014::MZ, pdg2014::Mt, 5 , 4);
      
      aMt6 = as5nf2as6nf(pdg2014::Mt, pdg2014::Mt, aMt5, /* nl= */5, 2);
      
      ak6 = run(aMt6, pdg2014::Mt, 100000, 6 );

      std::cout << " at5 = " << aMt5 << "  at6 = " << aMt6  << "  ak6 = " << ak6 << std::endl;

      std::cout <<  "3-loop decoupling " <<std::endl;
      aMt5 = run(pdg2014::asMZ, pdg2014::MZ, pdg2014::Mt, 5 , 4);
      
      aMt6 = as5nf2as6nf(pdg2014::Mt, pdg2014::Mt, aMt5, /* nl= */5, 3);
      
      ak6 = run(aMt6, pdg2014::Mt, 100000, 6 );

      std::cout << " at5 = " << aMt5 << "  at6 = " << aMt6  << "  ak6 = " << ak6 << std::endl;

      std::cout <<  "4-loop decoupling " <<std::endl;
      aMt5 = run(pdg2014::asMZ, pdg2014::MZ, pdg2014::Mt, 5 , 4);
      
      aMt6 = as5nf2as6nf(pdg2014::Mt, pdg2014::Mt, aMt5, /* nl= */5, 4);
      
      ak6 = run(aMt6, pdg2014::Mt, 100000, 6 );

      std::cout << " at5 = " << aMt5 << "  at6 = " << aMt6  << "  ak6 = " << ak6 << std::endl;

      std::cout << " aoi = " << aoi(50.) << std::endl;
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

