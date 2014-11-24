#include <iostream>
#include "mr.hpp"
// ROOT
// #include <TRandom1.h>
// #include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TApplication.h>
#include <TFitResult.h>

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Comare with:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0209084

      OSinput ACOVH80(0, 80.419, 91.188, 80, 174.3);
      OSinput ACOVH200(0, 80.419, 91.188, 200, 174.3);
      
      alphaMt  = 0.00779305;
      alphaS   = 0.1184;

      long double alphaMZ = 1./137.035999;

      HH<OS> dMH80  = HH<OS>(ACOVH80, ACOVH80.MMZ());
      HH<OS> dMH200 = HH<OS>(ACOVH200, ACOVH200.MMZ());
      
      // Compare with FIG.5
      std::cout << "1-loop \\alpha      Mh= " << ACOVH80.MH()  << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH80.x10(0,0) << std::endl;
      std::cout << "1-loop \\alpha      Mh= " << ACOVH200.MH() << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH200.x10(0,0) << std::endl;


      // Two-loop comaprison with Degrassi and Strumia and root-plot
      TApplication* rootapp = new TApplication("mH compare with ",&argc, argv);

      const Int_t n = 20;
      Double_t x[n];
      Double_t y[n];
      long double mustep  = 10; // GeV
      long double mustart = 80; // GeV

      AlphaS as;
            
      for (int mHi = 0; mHi < n; mHi++)
        {

          long double mH0 = 1.;

          long double MH = 0.01 + mH0 - mHi*mH0/double(n);

          OSinput DS2l(0, 80.384, 91.1876, MH , 173.10);

          HH<OS> dH  = HH<OS>(DS2l, DS2l.MMZ());          
          
          std::cout << " Mh= " << DS2l.MH() << ", [mH/MH -1] = " << dH.x20().real() << std::endl;     
          

          x[mHi] = MH;

          y[mHi] = 1./dH.x20().real();

        }

TGraph* gr = new TGraph(n,x,y);

// TF1  *mH2 = new TF1("mH2","1/x/x",0.01,1.);
// TF1  *mH4 = new TF1("mH4","1/x/x/x/x");

      // Trying to fit in mu=[80,180] range
      TFitResultPtr r = gr->Fit("pol2","S");
      // TFitResultPtr r = gr->Fit("mH2","S");
      
      std::cout << " b -: " << r->Value(0) << std::endl;
      gr->Draw();
      rootapp->Run();


      // // Two-loop comaprison with Degrassi and Strumia and root-plot
      // TApplication* rootapp = new TApplication("mH compare with ",&argc, argv);

      // const Int_t n = 10;
      // Double_t x[n];
      // Double_t y[n];
      // long double mHstep  = 10; // GeV
      // long double mHstart = 80; // GeV

      // AlphaS as;
            
      // for (int mHi = 0; mHi < n; mHi++)
      //   {
      //     OSinput DS2l(0, 80.384, 91.1876, mHstart + mHi*mHstep, 173.10);
      //     HH dH  = HH(DS2l, DS2l.MMt());          

      //     std::cout << "2-loop \\alpha^2      Mh= " << DS2l.MH() << ", [mH/MH -1] = " << DS2l.MMH()*pow(alphaMt/4./Pi,2)*dH.m20() << std::endl;     
          
      //     x[mHi] = DS2l.MH();
      //     y[mHi] = DS2l.MMH()*(
      //                          // alpha^2
      //                          // pow(alphaMt/4./Pi,2)*dH.m20().real() +
      //                          // alpha*lpha_S
      //                          alphaMt/4./Pi*as(DS2l.MMt())/4./Pi*dH.m11().real()
      //                          );
                               

      //   }

      // TGraph* gr = new TGraph(n,x,y);
      
      // gr->Draw();
      // rootapp->Run();
      
 
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

