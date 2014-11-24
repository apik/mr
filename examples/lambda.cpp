#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

long double alpha(long double mu)
{
  long double MZ  = 91.1876;
  long double aMZ = 1./127.944;
  return aMZ/(1+11./6./Pi*aMZ*log(mu/MZ));
}

int main (int argc, char *argv[])
{
  try
    {

      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS,alphaSMt;

      // Compare with:
      alphaMt   = 1./128.175;
      alphaSMt   = 0.1079;
      
     // Plot Yukawa top
      Plot4 plotLam("lam", "Higgs self-coupling", "mH", "delta_Lam", "a","a*a_S","a+a*a_S", "a+a*a_s+a^2");
      Plot3 plotDegr("lam2l", "Higgs self-coupling", "mH", "lam", "1+a+a*a_s","1+a+a*a_s+a^2","Degrassi");
      Plot4 plotDegrPart("lamPart", "Higgs self-coupling parts", "mH", "lam", "a*a_s","a^2","D:a*a_s","D:a^2");
      long double mHstep  = 10; // GeV
      long double mHstart = 110; // GeV
      
      OSinput BKKS(4.40, 80.384, 91.1876, mHstart, 173.10);     //Degrassi
      // OSinput BKKS(4.40, 80.385, 91.1876, mHstart, 173.5);   // Bezrukov
      // OSinput BKKS(4.40, 80.385, 91.1876, mHstart, 173.5);   // Kalmykov
      for (int mHi = 0; mHi < 4; mHi++)
        {
      
          BKKS.setMH(mHstart + mHi*mHstep);
          HH<OS> dH  = HH<OS>(BKKS, BKKS.MMt());          
          
          plotLam.add(BKKS.MH(),
                      alpha(BKKS.Mt())/4./Pi*dH.y10().real(),
                      alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dH.y11().real(),
                      alpha(BKKS.Mt())/4./Pi*dH.y10().real()+
                      alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dH.y11().real(),
                      alpha(BKKS.Mt())/4./Pi*dH.y10().real()+
                      alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dH.y11().real()+
                      pow(alpha(BKKS.Mt())/4./Pi,2)*dH.y20().real());

          plotDegr.add(BKKS.MH(),
                       (1+
                        alpha(BKKS.Mt())/4./Pi*dH.y10().real()+
                        alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dH.y11().real()// +
                        // pow(alpha(BKKS.Mt())/4./Pi,2)*dH.y20().real()
                        )*BKKS.MMH()*1.16637e-5/sqrt(2),
                       (1+
                        alpha(BKKS.Mt())/4./Pi*dH.y10().real()+
                        alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dH.y11().real()+
                        pow(alpha(BKKS.Mt())/4./Pi,2)*dH.y20().real())*BKKS.MMH()*1.16637e-5/sqrt(2),
                       (0.12711+0.00206*(BKKS.MH()-125.66)-0.00004*(BKKS.Mt()-173.10))
                       );
          
          plotDegrPart.add(BKKS.MH(),
                       alphaMt/4./Pi*alphaSMt/4./Pi*dH.y11().real()                  
                       *BKKS.MMH()*1.16637e-5/sqrt(2),
                       pow(alphaMt/4./Pi,2)*dH.y20().real()*BKKS.MMH()*1.16637e-5/sqrt(2),
                       (-23.89 + 0.12*(BKKS.MH()-125) - 0.64*(BKKS.Mt()-173))/pow(4*Pi,3)*alphaSMt,
                       (-9.45  - 0.11*(BKKS.MH()-125) - 0.21*(BKKS.Mt()-173))/pow(4*Pi,4)
                       );
          
        }
      
      
      

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

