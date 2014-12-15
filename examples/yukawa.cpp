#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS,alphaSMt;

      // Compare with:
      OSinput KVPhys(4.40, 80.385, 91.1876, 125.6, 173.5);
      

      alphaS   = 0.1184;

      long double alphaSMZ = 0.1184;
      // \mu = Mt
      alphaMt  = 0.00779305;
      alphaSMt   = 0.1184;

      
      // \mu = Mb
      long double alphaMb  = 0.00784257;
      long double alphaSMb   = 0.1905;

      

      AlphaS as;
      
      long double alphaMZ = 1./137.035999;


      // Yukawa top
      tt<OS> dMt(KVPhys, KVPhys.MMt());
      std::cout << "[ Top quark ]" << std::endl;
      std::cout << "Mh= " << KVPhys.MH()  << std::endl;
      std::cout << "as(MMt) = " << as(KVPhys.MMt()) << std::endl;          
      std::cout << "\t1-loop \\alpha         " << alphaMt/4./Pi*dMt.y10() << std::endl;
      std::cout << "\t1-loop \\alpha_S       " << alphaSMt/4./Pi*dMt.y01() << std::endl;
      std::cout << "\t2-loop \\alpha*\\alpha_S" << alphaMt/4./Pi*alphaSMt/4./Pi*dMt.y11() << std::endl;
      std::cout << "\t2-loop \\alpha^2       " << pow(alphaMt/4./Pi,2)*dMt.y20() << std::endl;


      // Yukawa bottom
      bb dMb  = bb(KVPhys, KVPhys.MMb());
      std::cout << "[ Bottom quark ]" << std::endl;
      std::cout << "\t1-loop \\alpha         " << alphaMb/4./Pi*dMb.y10() << std::endl;
      std::cout << "\t1-loop \\alpha_S       " << alphaSMb/4./Pi*dMb.y01() << std::endl;
      std::cout << "\t2-loop \\alpha*\\alpha_S" << alphaMb/4./Pi*alphaSMt/4./Pi*dMb.y11() << std::endl;
      std::cout << "\t2-loop \\alpha^2       " << pow(alphaMb/4./Pi,2)*dMb.y20() << std::endl;
      
      
      std::cout << "[ ratios ]" << std::endl;
      std::cout << "<10> yt/yb = " << dMb.y10()/dMt.y10() << " mt/mb = " << dMb.x10()/dMt.x10() << std::endl;
      std::cout << "<01> yt/yb = " << dMb.y01()/dMt.y01() << " mt/mb = " << dMb.x01()/dMt.x01() << std::endl;
      std::cout << "<11> yt/yb = " << dMb.y11()/dMt.y11() << " mt/mb = " << dMb.x11()/dMt.x11() << std::endl;
      std::cout << "<20> yt/yb = " << dMb.y20()/dMt.y20() << " mt/mb = " << dMb.x20()/dMt.x20() << std::endl;
      
      

      // Test Jegerlehner input
      // using 1-loop matching

      OSinput inFJ(0,80.385,91.1876,125.5,173.5);

      tt<OS> topFJ(inFJ, inFJ.MMt());

      std::cout << "\n\n \t Jegerlehner input:" << std::endl;
      
      std::complex<long double> ytFJ = inFJ.Mt()*(1 
                                                  + alphaMt/4./Pi*topFJ.x10() 
                                                  + alphaSMt/4./Pi*topFJ.x01()
                                                  + alphaMt/4./Pi*alphaSMt/4./Pi*topFJ.x11()
                                                  )*sqrt(2*sqrt(2)*1.16637e-5);

      std::cout << "yT(mt) = " << ytFJ << std::endl;


      // Plot Yukawa top
   Plot1 plotYt("yt", "Yukawa Top", "mH", "\\sigma_\\alpha*\\alpha_S", "a*a_S");
   long double mHstep  = 20; // GeV
   long double mHstart = 80; // GeV
   
   for (int mHi = 0; mHi < 13; mHi++)
     {
       OSinput DS2l(4.40, 80.385, 91.1876, mHstart + mHi*mHstep, 173.5);
       tt<OS> dtY(DS2l, DS2l.MMt());          
       
       plotYt.add(DS2l.MH(),alphaMt/4./Pi*alphaSMt/4./Pi*dtY.y11());
     }


   // Plot Yukawa bottom
   Plot1 plotYb("yb", "Yukawa Bottom", "mH", "\\sigma_\\alpha*\\alpha_S", "a*a_S");
   mHstep  = 10; // GeV
   mHstart = 1; // GeV
   
   for (int mHi = 0; mHi < 100; mHi++)
     {
       OSinput DS2l(4.40, 80.385, 91.1876, 125.6, 173.5);
       bb dbY  = bb(DS2l, mHstart + mHi*mHstep);          
       
       plotYb.add(mHstart + mHi*mHstep,alphaMb/4./Pi*alphaSMb/4./Pi*dbY.y11());
     }

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

