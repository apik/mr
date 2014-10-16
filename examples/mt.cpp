#include <iostream>
#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      long double xxx=1.;

      // Compare with:
      std::vector<OSinput> KV;
      KV.push_back(OSinput(0, xxx*80.385, xxx*91.1876, 124, 173.5));
      KV.push_back(OSinput(0, xxx*80.385, xxx*91.1876, 125, 173.5));
      KV.push_back(OSinput(0, xxx*80.385, xxx*91.1876, 126, 173.5));
      KV.push_back(OSinput(0, xxx*80.385, xxx*91.1876, 126, 173.5));
      KV.push_back(OSinput(0, xxx*80.385, xxx*91.1876, 126, 183.5));
      KV.push_back(OSinput(0, xxx*80.385, xxx*91.1876, 126, 193.5));
      
      alphaMt  = 0.00779305;
      alphaS   = 0.1184;

      AlphaS as;
      
      long double alphaMZ = 1./137.035999;


      for (std::vector<OSinput>::iterator it = KV.begin(); it != KV.end(); ++it)
        {
          tt dMt  = tt(*it, it->MMZ());
          std::cout << "Mh= " << it->MH()  << std::endl;
          std::cout << "as(MMt) = " << as(it->MMt()) << std::endl;          
          std::cout << "\t1-loop \\alpha         " << pow(xxx,2)*it->Mt()*alphaMt/4./Pi*dMt.m10() << std::endl;
          std::cout << "\t1-loop \\alpha_S       " << it->Mt()*as(it->MMt())/4./Pi*dMt.m01() << std::endl;
          std::cout << "\t2-loop \\alpha*\\alpha_S" << pow(xxx,2)*it->Mt()*alphaMt/4./Pi*// as(it->MMt())
            alphaS/4./Pi*dMt.m11() << std::endl;
          std::cout << "\t2-loop \\alpha^2  g.l. " << pow(xxx,4)*it->Mt()*pow(alphaMt/4./Pi,2)*dMt.mgl20() << std::endl;
          std::cout << "\t2-loop \\alpha^2       " << pow(xxx,4)*it->Mt()*pow(alphaMt/4./Pi,2)*dMt.m20() << std::endl;
          

          std::cout << " Full one-loop +QCD " << it->Mt()*alphaMt/4./Pi*dMt.m10() + it->Mt()*alphaS/4./Pi*dMt.m01() << std::endl;
          // dMt.test();
        }
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

