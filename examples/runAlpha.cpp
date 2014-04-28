#include <iostream>
#include <cmath>
#include "mr.hpp"
#include "CRunDec.h"

int main (int argc, char *argv[])
{
  try
    {
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

  
      MMt = pow(173.5,2);
      MMW = pow(80.385,2);
      MMZ = pow(91.1876,2);
      MMH = pow(125.66,2);
      
      alphaMt  = 0.00779305;
      alphaS = 0.1184;

      // Default loops=4, nf=5
      AlphaS as4l5;
      AlphaS as3l5(3,5);
      AlphaS as2l5(2,5);
      AlphaS as1l5(1,5);

      // CRunDec for comparison nf=5
      CRunDec* pObjnf5 = new CRunDec(5);      

      std::cout.precision(20);
      std::cout << std::fixed;

      std::cout << "loops = 4, nf = 5, \\alpha_s(MMt) = " << as4l5(MMt) << std::endl;

      long double mu20 = MMZ;
      long double step = (MMt - MMZ)/10;
      for(int i = 1; i < 10; i++)
        {
          long double mu2 = mu20 + i*step;

          std::cout << "Difference for \\mu=" << sqrt(mu2) << std::endl 
                    << "\tin 4-loop: " << as4l5(mu2) - pObjnf5 -> AlphasExact(alphaS,sqrt(MMZ),sqrt(mu2),4) << " "
                    << 100*fabs((as4l5(mu2) - pObjnf5 -> AlphasExact(alphaS,sqrt(MMZ),sqrt(mu2),4))/as4l5(mu2)) << " % " << std::endl

                    << "\tin 3-loop: " << as3l5(mu2) - pObjnf5 -> AlphasExact(alphaS,sqrt(MMZ),sqrt(mu2),3) << " "
                    << 100*fabs((as3l5(mu2) - pObjnf5 -> AlphasExact(alphaS,sqrt(MMZ),sqrt(mu2),3))/as3l5(mu2)) << " % "<< std::endl

                    << "\tin 2-loop: " << as2l5(mu2) - pObjnf5 -> AlphasExact(alphaS,sqrt(MMZ),sqrt(mu2),2) << " "
                    << 100*fabs((as2l5(mu2) - pObjnf5 -> AlphasExact(alphaS,sqrt(MMZ),sqrt(mu2),2))/as2l5(mu2)) << " % "<< std::endl

                    << "\tin 1-loop: " << as1l5(mu2) - pObjnf5 -> AlphasExact(alphaS,sqrt(MMZ),sqrt(mu2),1) << " "
                    << 100*fabs((as1l5(mu2) - pObjnf5 -> AlphasExact(alphaS,sqrt(MMZ),sqrt(mu2),1))/as1l5(mu2)) << " % "<< std::endl << std::endl;

        }

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

