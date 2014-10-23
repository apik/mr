// #include <CRunDec.h>
#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"



long double alpha(long double mu)
{
  long double MZ  = 91.1876;
  long double aMZ = 1./127.944;
  return aMZ/(1-11./6./Pi*aMZ*log(mu/MZ));
}

int main (int argc, char *argv[])
{
  try
    {

      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS,alphaSMt;

      // Scale inv test
      long double a = 1.;
      // Compare with:
      alphaMt   = 1./128.175;
      alphaSMt   = 0.1079;
      
      std::vector<OSinput> sv;

      sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125.6, a*173.5)); // A
      // sv.push_back(OSinput(a*4.4, a*80.399, a*91.1876, a*126, a*173.9)); // B
      // sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125, a*172.9)); // C

      // sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125, a*274.1)); // C
      // sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125, a*374.1)); // C
      // sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125, a*474.1)); // C
      


      for (std::vector<OSinput>::iterator it = sv.begin(); it != sv.end(); ++it)
        {
          tt dMt  = tt(*it, it->MMZ());
          HH dMH  = HH(*it, it->MMZ());
          dr ddr  = dr(*it, it->MMZ());
          std::cout << "sigmata   [1,0] = " << dMt.m10() << std::endl;
          std::cout << "deltayta  [1,0]= " << dMt.my10() << std::endl;
          std::cout << "sigmataaS [1,1]= " << dMt.m11() << std::endl;
          std::cout << "deltaytaaS[1,1] = " << dMt.my11() << std::endl;
          
          std::cout << "\n\n        dr[1,0] = " << ddr.dr10() << std::endl;
          std::cout << "\n\n        dr[1,1] = " << ddr.dr11() << std::endl;

          std::cout << "sigmaHa   [1,0] = " << dMH.m10() << std::endl;
          std::cout << "deltayHa  [1,0]= " << dMH.my10() << std::endl;
          std::cout << "sigmaHaaS [1,1]= " << dMH.m11() << std::endl;
          std::cout << "deltayHaaS[1,1] = " << dMH.my11() << std::endl;



        }
      
      

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

