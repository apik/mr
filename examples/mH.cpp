#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

  
      MMt = pow(173.5,2);
      MMW = pow(80.385,2);
      MMZ = pow(91.1876,2);
      MMH = pow(125.66,2);
      
      alphaMt  = 0.00779305;
      alphaS = 0.1184;

      HH h2l = HH(MMt, MMH, MMW, MMZ, MMt);

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

