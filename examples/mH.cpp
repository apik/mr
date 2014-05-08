#include <iostream>
#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Comare with:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0209084

      SMinput ACOVH80(80.419, 91.188, 80, 174.3);
      SMinput ACOVH200(80.419, 91.188, 200, 174.3);
      
      alphaMt  = 0.00779305;
      alphaS   = 0.1184;

      long double alphaMZ = 1./137.035999;

      HH dMH80  = HH(ACOVH80, ACOVH80.MMZ());
      HH dMH200 = HH(ACOVH200, ACOVH200.MMZ());
      
      // Compare with FIG.5
      std::cout << "1-loop \\alpha      Mh= " << ACOVH80.mh()  << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH80.m10(0,0) << std::endl;
      std::cout << "1-loop \\alpha      Mh= " << ACOVH200.mh() << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH200.m10(0,0) << std::endl;
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

