#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Compare with:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0209084

      // OSinput ACOVH80(0, 80.419, 91.188, 80, 174.3);
      // OSinput ACOVH200(0, 80.419, 91.188, 200, 174.3);
      
      // alphaMt  = 0.00779305;
      // alphaS   = 0.1185;

      // // long double alphaMZ = 1./137.035999;

      long double alphaMZ = 1./128.992;
      // HH<OS> dMH80  = HH<OS>(ACOVH80, ACOVH80.MMZ());
      // HH<OS> dMH200 = HH<OS>(ACOVH200, ACOVH200.MMZ());
      
      // // Compare with FIG.5
      // std::cout << "1-loop \\alpha      Mh= " << ACOVH80.MH()  << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH80.m10(0,0) << std::endl;
      // std::cout << "1-loop \\alpha      Mh= " << ACOVH200.MH() << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH200.m10(0,0) << std::endl;

      BetaVEV bv(3);

      double vi[8] = {1.,2.,3.,4.,5.,6.,7.,8.}; 
      state_type ain(vi, vi + 8);

      
      
      std::cout << "VEV = " << BetaVEV::gamv(ain,3) << std::endl;
      std::cout << "bms = " << BetaMu2::bmu2(ain,3) << std::endl;

      
      long double Mt = 173.10;
      // MSinput SPM = MSinput::fromConsts( 0.127,
      //                                    367.241, //247.0
      //                                    0, 0.936, 0.648, 0.358);

      // std::cout << "Mb= " << SPM.mb() << std::endl; 
      // std::cout << "MW= " << SPM.mW() << std::endl; 
      // std::cout << "MZ= " << SPM.mZ() << std::endl; 
      // std::cout << "MH= " << SPM.mH() << std::endl; 
      // std::cout << "Mt= " << SPM.mt() << std::endl; 
      // std::cout << "1/al= " << 1./SPM.alpha() << std::endl; 

      
      // HH<MS> mHH(SPM, pow(Mt,2));
      // std::cout << "1-loop \\alpha      Mh= " << SPM.mH() << ", [mH/MH -1] = " << SPM.mH()*sqrt(1 + SPM.alpha()/4./Pi*mHH.x10()) << std::endl;


      // OSinput SMH(0, 80.419, 91.188, 125.6, 172.3);


      // HH<OS> dHHH(SMH, SMH.MMt());

      // std::cout << "1-loop \\alpha   [mH(mu)] = "  << sqrt(SMH.MMH()*(1+alphaMZ/4./Pi*dHHH.x10())) << std::endl;

      // HH<OS> dHHH_Planck(SMH, pow(1.2209,2) * pow(10.,2*19));
      // double  alphaMpl = 1./104.82;
      // double  alphaSMpl = 0.0188966;
      // std::cout << "1-loop \\alpha   [mH(mPl)] = "  << sqrt(SMH.MMH()*(1
      //                                                                  + alphaMpl/4./Pi*dHHH.x10()
      //                                                                  + alphaMpl/4./Pi*alphaSMpl/4./Pi*dHHH.x11())) << std::endl;

      // std::cout << "1-loop \\alpha   [mH/MH -1] = "  << SMH.MMH()*(1+alphaMZ/4./Pi*dHHH.x10()) << std::endl;
      // std::cout << "1-loop \\alpha               Mh= " << BKKS2l200.MH() << ", [mH/MH -1] = "  << alphaMZ/4./Pi*dEW200.m10() << std::endl;

      // std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << BKKS2l80.MH()  << ", [mH/MH -1] = "  << alphaMZ/4./Pi*alphaS/4./Pi*dEW80.m11() << std::endl;
      // std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << BKKS2l200.MH() << ", [mH/MH -1] = "  << alphaMZ/4./Pi*alphaS/4./Pi*dEW200.m11() << std::endl;

      
      // std::cout << "2-loop \\alpha      Mh= " << SPM.mH() << ", [mH/MH -1] = " << pow(SPM.alpha()/4./Pi,1)*mHH.m10() << std::endl;

      // // 
      // // Compare with BKKS arXiv:1205.2893           
      // // 
      // std::cout << "BKKS arXiv:1205.2893" << std::endl;

      // OSinput BKKS2l80 (0, 80.419, 91.188, 80, 174.3);
      // OSinput BKKS2l200(0, 80.419, 91.188, 200, 174.3);


      // HH dEW80  = HH(BKKS2l80, BKKS2l80.MMt());
      // HH dEW200 = HH(BKKS2l200, BKKS2l200.MMt());
      
      // std::cout << "1-loop \\alpha               Mh= " << BKKS2l80.MH()  << ", [mH/MH -1] = "  << alphaMZ/4./Pi*dEW80.m10() << std::endl;
      // std::cout << "1-loop \\alpha               Mh= " << BKKS2l200.MH() << ", [mH/MH -1] = "  << alphaMZ/4./Pi*dEW200.m10() << std::endl;

      // std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << BKKS2l80.MH()  << ", [mH/MH -1] = "  << alphaMZ/4./Pi*alphaS/4./Pi*dEW80.m11() << std::endl;
      // std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << BKKS2l200.MH() << ", [mH/MH -1] = "  << alphaMZ/4./Pi*alphaS/4./Pi*dEW200.m11() << std::endl;


      // Plot3 plot("mH", "\\alpha*\\alpha_S conversion between \\bar{MS} and  On-Shell Higgs mass", "mH", "m/M-1", "1-loop ", "2-loop EM*QCD", "2-loop EM^2");
      // AlphaS as;
      // long double mHstep  = 10; // GeV
      // long double mHstart = 80; // GeV

      // for (int mHi = 0; mHi < 13; mHi++)
      //   {
      //     OSinput DS2l(0, 80.384, 91.1876, mHstart + mHi*mHstep, 173.10);
      //     HH dH  = HH(DS2l, DS2l.MMt());          
          
      //     plot.add(DS2l.MH(), 
      //              alphaMt/4./Pi*dH.m10().real(), 
      //              alphaMt/4./Pi*as(DS2l.MMt())/4./Pi*dH.m11().real(),
      //              pow(alphaMt/4./Pi,2)*dH.m20().real());
      //   }

      // // 
      // // Two-loop comaprison with Degrassi and Strumia
      // // 
      // OSinput DS2l80(0, 80.384, 91.1876, 80, 173.10);
      // OSinput DS2l200(0, 80.384, 91.1876, 245, 173.10);


      // HH dEWEW80  = HH(DS2l80, DS2l80.MMt());
      // HH dEWEW200 = HH(DS2l200, DS2l200.MMt());
      
      // std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << DS2l80.MH() 
      //           << ", [mH/MH -1] = " << DS2l80.MMH()*alphaMt/4./Pi*as(DS2l80.MMt())/4./Pi*dEWEW80.m11() << std::endl;
      // std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << DS2l200.MH() 
      //           << ", [mH/MH -1] = " << DS2l200.MMH()*alphaMt/4./Pi*as(DS2l200.MMt())/4./Pi*dEWEW200.m11() << std::endl;

      // std::cout << "2-loop \\alpha^2      Mh= " << DS2l80.MH() 
      //           << ", [mH/MH -1] = " << DS2l80.MMH()*pow(alphaMt/4./Pi,2)*dEWEW80.m20() << std::endl;
      // std::cout << "2-loop \\alpha^2      Mh= " << DS2l200.MH() 
      //           << ", [mH/MH -1] = " << DS2l200.MMH()*pow(alphaMt/4./Pi,2)*dEWEW200.m20() << std::endl;

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

