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
      long double a = 40.;
      // Compare with:
      alphaMt   = 1./128.175;
      alphaSMt   = 0.1079;
      
      std::vector<SMinput> sv;

      sv.push_back(SMinput(a*0, a*80.399, a*91.1876, a*125, a*173.1)); // A
      sv.push_back(SMinput(a*4.4, a*80.399, a*91.1876, a*126, a*173.1)); // B
      sv.push_back(SMinput(a*0, a*80.399, a*91.1876, a*125, a*174.1)); // C

      sv.push_back(SMinput(a*0, a*80.399, a*91.1876, a*125, a*274.1)); // C
      sv.push_back(SMinput(a*0, a*80.399, a*91.1876, a*125, a*374.1)); // C
      sv.push_back(SMinput(a*0, a*80.399, a*91.1876, a*125, a*474.1)); // C
      


      for (std::vector<SMinput>::iterator it = sv.begin(); it != sv.end(); ++it)
        {
          tt dMt  = tt(*it, it->MMt());
         
          std::cout << "dytalas = " << dMt.my11() << std::endl;

        }
      
      

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

