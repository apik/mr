#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Comare with:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0105304
      // fermionic    part (nH = 1, nL = 2) : arXiv:hep-ph/0212319           
      SMinput JKVH80(80.419, 91.188, 80, 174.3);
      SMinput JKVH200(80.419, 91.188, 200, 174.3);
         
      alphaMt  = 0.00779305;
      alphaS = 0.1184;

      long double alphaMZ = 1./137.035999;

      WW W80  = WW(JKVH80, JKVH80.MMW());
      WW W200 = WW(JKVH200, JKVH200.MMW());

      std::cout << " SM input MMZ: " << JKVH200.MMt() << std::endl;
      // 
      // Bosonic part
      // 
      std::cout << alphaMZ/4./Pi*W80.m10(0,0) << std::endl;
      std::cout << alphaMZ/4./Pi*W200.m10(0,0) << std::endl;

      // 
      // Full SM
      // 
      std::cout << alphaMZ/4./Pi*W80.m10(2,1) << std::endl;
      std::cout << alphaMZ/4./Pi*W200.m10(2,1) << std::endl;


    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

