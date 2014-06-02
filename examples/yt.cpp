// #include <CRunDec.h>
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
      Plot5 plotYt("ytQCD", "Yukawa top", "mH", "delta_Yt", "as^2","as^3","as^3+a", "as^3+a+a*as","as^3+a+a*as+a^2" );
      // Plot3 plotYtDegr("yt_degr", "Higgs self-coupling", "mH", "lam", "1+a+a*a_s","1+a+a*a_s+a^2","Degrassi");
      Plot4 plotYtDegrPart("ytPart", "Yukawa top parts", "mH", "lam", "a*a_s","a^2","D:a*a_s","D:a^2");
      long double mHstep  = 10; // GeV
      long double mHstart = 110; // GeV
      
      SMinput BKKS(4.40, 80.384, 91.1876, mHstart, 173.10);
      // SMinput BKKS(4.40, 80.385, 91.1876, mHstart, 173.5);   // Bezrukov
      // SMinput BKKS(4.40, 80.385, 91.1876, mHstart, 173.5);   // Kalmykov
      for (int mHi = 0; mHi < 4; mHi++)
        {
      
          BKKS.setMH(mHstart + mHi*mHstep);
          tt dt  = tt(BKKS, BKKS.MMt());          
          
          plotYt.add(BKKS.MH(),
                     pow(alphaSMt/4./Pi,1)*dt.my01().real()+
                     pow(alphaSMt/4./Pi,2)*dt.m02().real(),
                     
                     pow(alphaSMt/4./Pi,1)*dt.my01().real()+
                     pow(alphaSMt/4./Pi,2)*dt.m02().real()+
                     pow(alphaSMt/4./Pi,3)*dt.m03().real(),

                     pow(alphaSMt/4./Pi,1)*dt.m01().real()+
                     pow(alphaSMt/4./Pi,2)*dt.m02().real()+
                     pow(alphaSMt/4./Pi,3)*dt.m03().real()+
                     alpha(BKKS.Mt())/4./Pi*dt.my10().real(),

                     pow(alphaSMt/4./Pi,1)*dt.m01().real()+
                     pow(alphaSMt/4./Pi,2)*dt.m02().real()+
                     pow(alphaSMt/4./Pi,3)*dt.m03().real()+
                     alpha(BKKS.Mt())/4./Pi*dt.my10().real()+
                     alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dt.my11().real(),

                     pow(alphaSMt/4./Pi,1)*dt.m01().real()+
                     pow(alphaSMt/4./Pi,2)*dt.m02().real()+
                     pow(alphaSMt/4./Pi,3)*dt.m03().real()+
                     alpha(BKKS.Mt())/4./Pi*dt.my10().real()+
                     alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dt.my11().real()+
                     pow(alpha(BKKS.Mt())/4./Pi,2)*dt.my20().real());

          // plotDegr.add(BKKS.MH(),
          //              (1+
          //               alpha(BKKS.Mt())/4./Pi*dH.lam10().real()+
          //               alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dH.lam11().real()// +
          //               // pow(alpha(BKKS.Mt())/4./Pi,2)*dH.lam20().real()
          //               )*BKKS.MMH()*1.16637e-5/sqrt(2),
          //              (1+
          //               alpha(BKKS.Mt())/4./Pi*dH.lam10().real()+
          //               alpha(BKKS.Mt())/4./Pi*alphaSMt/4./Pi*dH.lam11().real()+
          //               pow(alpha(BKKS.Mt())/4./Pi,2)*dH.lam20().real())*BKKS.MMH()*1.16637e-5/sqrt(2),
          //              (0.12711+0.00206*(BKKS.MH()-125.66)-0.00004*(BKKS.Mt()-173.10))
          //              );
          
          plotYtDegrPart.add(BKKS.MH(),
                             alphaMt/4./Pi*alphaSMt/4./Pi*dt.my11().real()*BKKS.Mt()*sqrt(1.16637e-5*sqrt(8)),
                             pow(alphaMt/4./Pi,2)*dt.my20().real()*BKKS.Mt()*sqrt(1.16637e-5*sqrt(8)),
                             (-7.53 + 0.09*(BKKS.MH()-125) - 0.23*(BKKS.Mt()-173))/pow(4*Pi,3)*alphaSMt,
                             (5.22  - 0.01*(BKKS.MH()-125) + 0.15*(BKKS.Mt()-173))/pow(4*Pi,4)
                             );
          
          // std::cout << "Rundec: " << 4.*alphaSMt/4./Pi*CRunDec::fMsFromOs1(BKKS.Mt(), BKKS.Mt()) << std::endl; 
          std::cout << "      : " << alphaSMt/4./Pi*dt.my01().real() << "  A= " << Tsil::A(BKKS.MMt(),BKKS.MMt())/BKKS.MMt() << std::endl; 

          
          
        }
      
      

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

