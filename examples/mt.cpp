#include <iostream>
#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Compare with:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0105304
      // fermionic    part (nH = 1, nL = 2) : arXiv:hep-ph/0212319           
      SMinput* KV[3];
      KV[0] = new SMinput(80.385, 91.1876, 124, 173.5);
      KV[1] = new SMinput(80.385, 91.1876, 125, 173.5);
      KV[2] = new SMinput(80.385, 91.1876, 126, 173.5);
      
      alphaMt  = 0.00779305;
      alphaS   = 0.1184;

      AlphaS as;
      
      long double alphaMZ = 1./137.035999;


      
      for (int i = 0; i < 3; i++)
        {
          tt dMt  = tt(*KV[i], KV[i]->MMt());
          std::cout << "Mh= " << KV[i]->mh()  << std::endl;
          std::cout << "as(MMt) = " << as(KV[i]->MMt()) << std::endl;          
          std::cout << "\t1-loop \\alpha         " << KV[i]->mt()*alphaMt/4./Pi*dMt.m10() << std::endl;
          std::cout << "\t1-loop \\alpha_S       " << KV[i]->mt()*as(KV[i]->MMt())/4./Pi*dMt.m01() << std::endl;
          std::cout << "\t2-loop \\alpha*\\alpha_S" << KV[i]->mt()*alphaMt/4./Pi*as(KV[i]->MMt())/4./Pi*dMt.m11() << std::endl;
          std::cout << "\t2-loop \\alpha^2       " << KV[i]->mt()*pow(alphaMt/4./Pi,2)*dMt.m20() << std::endl;
        }
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

