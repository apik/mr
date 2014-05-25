#include <iostream>
#include "mr.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      // fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Compare with:
      SMinput* KV[3];
      KV[0] = new SMinput(4.40204, 80.385, 91.1876, 124, 173.5);
      KV[1] = new SMinput(4.40204, 80.385, 91.1876, 125, 173.5);
      KV[2] = new SMinput(4.40204, 80.385, 91.1876, 126, 173.5);
      
      alphaMt  = 0.00779305;
      alphaS   = 0.1184;

      AlphaS as;
      
      long double alpha0 = 1./137.035999;

      long double asMb =  0.1905;

      
      for (int i = 0; i < 3; i++)
        {
          bb dMb  = bb(*KV[i], KV[i]->MMb());
          std::cout << "Mh= " << KV[i]->MH()  << std::endl;
          std::cout << "as(MMt) = " << as(KV[i]->MMt()) << std::endl;          
          std::cout << "\t1-loop \\alpha         " << // KV[i]->Mb()*
            alpha0/4./Pi*dMb.m10() << std::endl;
          std::cout << "\t1-loop \\alpha_S       " << // KV[i]->Mb()*
            asMb/4./Pi*dMb.m01() << std::endl;
          std::cout << "\t2-loop \\alpha*\\alpha_S" << // KV[i]->Mb()*
            alpha0/4./Pi*asMb/4./Pi*dMb.m11() << std::endl;
          std::cout << "\t2-loop \\alpha^2       " << // KV[i]->Mb()*
            pow(alpha0/4./Pi,2)*dMb.m20() << std::endl;
        }
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

