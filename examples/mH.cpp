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

      // Comare with:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0209084

      SMinput ACOVH80(80.419, 91.188, 80, 174.3);
      SMinput ACOVH200(80.419, 91.188, 200, 174.3);
      
      alphaMt  = 0.00779305;
      alphaS   = 0.1185;

      // long double alphaMZ = 1./137.035999;

      long double alphaMZ = 1./128.992;
      HH dMH80  = HH(ACOVH80, ACOVH80.MMZ());
      HH dMH200 = HH(ACOVH200, ACOVH200.MMZ());
      
      // Compare with FIG.5
      std::cout << "1-loop \\alpha      Mh= " << ACOVH80.MH()  << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH80.m10(0,0) << std::endl;
      std::cout << "1-loop \\alpha      Mh= " << ACOVH200.MH() << ", [mH/MH -1] = " << alphaMZ/4./Pi*dMH200.m10(0,0) << std::endl;
      
      // 
      // Compare with BKKS arXiv:1205.2893           
      // 
      std::cout << "BKKS arXiv:1205.2893" << std::endl;

      SMinput BKKS2l80 (80.419, 91.188, 80, 174.3);
      SMinput BKKS2l200(80.419, 91.188, 200, 174.3);


      HH dEW80  = HH(BKKS2l80, BKKS2l80.MMt());
      HH dEW200 = HH(BKKS2l200, BKKS2l200.MMt());
      
      std::cout << "1-loop \\alpha               Mh= " << BKKS2l80.MH()  << ", [mH/MH -1] = "  << alphaMZ/4./Pi*dEW80.m10() << std::endl;
      std::cout << "1-loop \\alpha               Mh= " << BKKS2l200.MH() << ", [mH/MH -1] = "  << alphaMZ/4./Pi*dEW200.m10() << std::endl;

      std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << BKKS2l80.MH()  << ", [mH/MH -1] = "  << alphaMZ/4./Pi*alphaS/4./Pi*dEW80.m11() << std::endl;
      std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << BKKS2l200.MH() << ", [mH/MH -1] = "  << alphaMZ/4./Pi*alphaS/4./Pi*dEW200.m11() << std::endl;


      Plot3 plot("mH", "\\alpha*\\alpha_S conversion between \\bar{MS} and  On-Shell Higgs mass", "mH", "m/M-1", "1-loop ", "2-loop EM*QCD", "2-loop EM^2");
      AlphaS as;
      long double mHstep  = 10; // GeV
      long double mHstart = 80; // GeV

      for (int mHi = 0; mHi < 13; mHi++)
        {
          SMinput DS2l(80.384, 91.1876, mHstart + mHi*mHstep, 173.10);
          HH dH  = HH(DS2l, DS2l.MMt());          
          
          plot.add(DS2l.MH(), 
                   alphaMt/4./Pi*dH.m10().real(), 
                   alphaMt/4./Pi*as(DS2l.MMt())/4./Pi*dH.m11().real(),
                   pow(alphaMt/4./Pi,2)*dH.m20().real());
        }

      // 
      // Two-loop comaprison with Degrassi and Strumia
      // 
      SMinput DS2l80(80.384, 91.1876, 80, 173.10);
      SMinput DS2l200(80.384, 91.1876, 245, 173.10);


      HH dEWEW80  = HH(DS2l80, DS2l80.MMt());
      HH dEWEW200 = HH(DS2l200, DS2l200.MMt());
      
      std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << DS2l80.MH() 
                << ", [mH/MH -1] = " << DS2l80.MMH()*alphaMt/4./Pi*as(DS2l80.MMt())/4./Pi*dEWEW80.m11() << std::endl;
      std::cout << "2-loop \\alpha*\\alpha_S      Mh= " << DS2l200.MH() 
                << ", [mH/MH -1] = " << DS2l200.MMH()*alphaMt/4./Pi*as(DS2l200.MMt())/4./Pi*dEWEW200.m11() << std::endl;

      std::cout << "2-loop \\alpha^2      Mh= " << DS2l80.MH() 
                << ", [mH/MH -1] = " << DS2l80.MMH()*pow(alphaMt/4./Pi,2)*dEWEW80.m20() << std::endl;
      std::cout << "2-loop \\alpha^2      Mh= " << DS2l200.MH() 
                << ", [mH/MH -1] = " << DS2l200.MMH()*pow(alphaMt/4./Pi,2)*dEWEW200.m20() << std::endl;

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

