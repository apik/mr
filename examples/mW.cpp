#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Comare with JKVH:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0105304
      // fermionic    part (nH = 1, nL = 2) : arXiv:hep-ph/0212319           
      SMinput JKVH80(80.419, 91.188, 80, 174.3);
      SMinput JKVH200(80.419, 91.188, 200, 174.3);
         
      alphaMt  = 0.00779305;
      alphaS = 0.1184;

      long double alphaMZ = 1./137.035999;

      WW W80  = WW(JKVH80, JKVH80.MMW());
      WW W200 = WW(JKVH200, JKVH200.MMW());

      std::cout << "\\mu = MW" << std::endl;
      // 
      // Bosonic part 1-loop
      // 
      std::cout << "No fermions: 1-loop \\alpha      Mh= " << JKVH80.mh()  << ", [mW/MW -1] = " 
                << alphaMZ/4./Pi*W80.m10(0,0) << std::endl;
      std::cout << "No fermions: 1-loop \\alpha      Mh= " << JKVH80.mh()  << ", [mW/MW -1] = " 
                << alphaMZ/4./Pi*W200.m10(0,0) << std::endl;

      // 
      // Full SM
      // 
      std::cout << "SM 1-loop \\alpha      Mh= " << JKVH80.mh()  << ", [mW/MW -1] = " 
                << alphaMZ/4./Pi*W80.m10(2,1) << std::endl;
      std::cout << "SM 1-loop \\alpha      Mh= " << JKVH80.mh()  << ", [mW/MW -1] = " 
                << alphaMZ/4./Pi*W200.m10(2,1) << std::endl;


      // In ACOV \mu^2 = MMZ
      WW W_ACOV_80  = WW(JKVH80, JKVH80.MMZ());
      WW W_ACOV_200 = WW(JKVH200, JKVH200.MMZ());
      
      // 
      // Bosonic part 2-loop, compare with ACOV: arXiv:hep-ph/0209084,
      // FIG. 6
      // 
      std::cout << "\\mu = MZ" << std::endl;
      std::cout << "2-loop \\alpha^2      Mh= " << JKVH80.mh()  << ", [mW/MW -1] = " 
                << pow(alphaMZ/4./Pi,2)*W_ACOV_80.m20(0,0) << std::endl;
      std::cout << "2-loop \\alpha^2      Mh= " << JKVH200.mh()  << ", [mW/MW -1] = " 
                << pow(alphaMZ/4./Pi,2)*W_ACOV_200.m20(0,0) << std::endl;


    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

