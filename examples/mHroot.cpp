#include <iostream>
#include "mr.hpp"
// ROOT
// #include <TRandom1.h>
// #include <TH1F.h>
#include <TGraph.h>
#include <TApplication.h>

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Comare with:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0209084

      SMinput ACOVH80(0, 80.419, 91.188, 80, 174.3);
      SMinput ACOVH200(0, 80.419, 91.188, 200, 174.3);
      
      alphaMt  = 0.00779305;
      alphaS   = 0.1184;

      long double alphaMZ = 1./137.035999;

      HH dMH80  = HH(ACOVH80, ACOVH80.MMZ());
      HH dMH200 = HH(ACOVH200, ACOVH200.MMZ());
      
      // Compare with FIG.5
      std::cout << "1-loop \\alpha      Mh= " << ACOVH80.MH()  << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH80.m10(0,0) << std::endl;
      std::cout << "1-loop \\alpha      Mh= " << ACOVH200.MH() << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH200.m10(0,0) << std::endl;


      // Two-loop comaprison with Degrassi and Strumia and root-plot
      TApplication* rootapp = new TApplication("mH compare with ",&argc, argv);

      const Int_t n = 10;
      Double_t x[n];
      Double_t y[n];
      long double mHstep  = 10; // GeV
      long double mHstart = 80; // GeV

      AlphaS as;
            
      for (int mHi = 0; mHi < n; mHi++)
        {
          SMinput DS2l(0, 80.384, 91.1876, mHstart + mHi*mHstep, 173.10);
          HH dH  = HH(DS2l, DS2l.MMt());          

          std::cout << "2-loop \\alpha^2      Mh= " << DS2l.MH() << ", [mH/MH -1] = " << DS2l.MMH()*pow(alphaMt/4./Pi,2)*dH.m20() << std::endl;     
          
          x[mHi] = DS2l.MH();
          y[mHi] = DS2l.MMH()*(
                               // alpha^2
                               // pow(alphaMt/4./Pi,2)*dH.m20().real() +
                               // alpha*lpha_S
                               alphaMt/4./Pi*as(DS2l.MMt())/4./Pi*dH.m11().real()
                               );
                               

        }

      TGraph* gr = new TGraph(n,x,y);
      
      gr->Draw();
      rootapp->Run();
      
 
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

