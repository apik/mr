#include "mr.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      // Comare with JKVH:
      // bosonic only part (nH = nL = 0)    : arXiv:hep-ph/0105304
      // fermionic    part (nH = 1, nL = 2) : arXiv:hep-ph/0212319           
      SMinput JKVH80(0, 80.419, 91.188, 150, 174.3);
      SMinput JKVH200(0, 80.419, 91.188, 300, 174.3);
         
      alphaMt  = 0.00779305;
      alphaS = 0.1184;

      long double alpha0 = 1./137.035999;

      WW W80  = WW(JKVH80, JKVH80.MMW());
      WW W200 = WW(JKVH200, JKVH200.MMW());

      std::cout << "\\mu = MW" << std::endl;
      // 
      // Bosonic part 1-loop
      // 
      std::cout << "No fermions: 1-loop \\alpha      Mh= " << JKVH80.MH()  << ", [mW/MW -1] = " 
                << alpha0/4./Pi*W80.m10(0,0) << std::endl;
      std::cout << "No fermions: 1-loop \\alpha      Mh= " << JKVH200.MH()  << ", [mW/MW -1] = " 
                << alpha0/4./Pi*W200.m10(0,0) << std::endl;

      // 
      // Full SM
      // 
      std::cout << "SM 1-loop \\alpha      Mh= " << JKVH80.MH()  << ", [mW/MW -1] = " 
                << alpha0/4./Pi*W80.m10(2,1) << std::endl;
      std::cout << "SM 1-loop \\alpha      Mh= " << JKVH200.MH()  << ", [mW/MW -1] = " 
                << alpha0/4./Pi*W200.m10(2,1) << std::endl;
      std::cout << "SM 2-loop \\alpha^2 nf=0     Mh= " << JKVH80.MH()  << ", [mW/MW -1] = " 
                << pow(alpha0/4./Pi,2)*W80.m20(0,0) << std::endl;
      std::cout << "SM 2-loop \\alpha^2 nf=0     Mh= " << JKVH200.MH()  << ", [mW/MW -1] = " 
                << pow(alpha0/4./Pi,2)*W200.m20(0,0) << std::endl;

      std::cout << " 1 = (1 + X1 + X2)*(1 + Y1 + Y2)" << std::endl;
      
      std::cout << "LHS      Mh= " << JKVH80.MH()  << ", [mW/MW -1] = " 
                << 1 + alpha0/4./Pi*W80.m10() + 0*pow(alpha0/4./Pi,2)*W80.m20() << std::endl;

      std::cout << "LHS      Mh= " << JKVH200.MH()  << ", [mW/MW -1] = " 
                << 1 + alpha0/4./Pi*W200.m10() + 0*pow(alpha0/4./Pi,2)*W200.m20() << std::endl;

      std::cout << "LHS 1l     Mh= " << JKVH80.MH()  << ", [mW/MW -1] = " 
                << 0.9652791210924665318293049331952014859085794548969162728818377550731344617939375411461317952286357517 + pow(alpha0/4./Pi,2)*W80.m20() << std::endl;

      std::cout << "LHS 1l     Mh= " << JKVH200.MH()  << ", [mW/MW -1] = " 
                << 1.044246369785103344770588814451876892018947264237810660355287674109945054514877377628523175447663686 + pow(alpha0/4./Pi,2)*W200.m20() << std::endl;



      std::cout << "SM 1-loop \\alpha      Mh= " << JKVH80.MH()  << ", [MW/mW -1] = " 
                << alpha0/4./Pi*W80.m10(2,1)/(1+alpha0/4./Pi*W80.m10(2,1)) << std::endl;
      std::cout << "SM 1-loop \\alpha      Mh= " << JKVH200.MH()  << ", [MW/mW -1] = " 
                << alpha0/4./Pi*W200.m10(2,1)/(1+alpha0/4./Pi*W200.m10(2,1)) << std::endl;




      // In ACOV \mu^2 = MMZ
      WW W_ACOV_80  = WW(JKVH80, JKVH80.MMZ());
      WW W_ACOV_200 = WW(JKVH200, JKVH200.MMZ());
      
      // 
      // Bosonic part 2-loop, compare with ACOV: arXiv:hep-ph/0209084,
      // FIG. 6
      // 
      std::cout << "\\mu = MZ" << std::endl;
      std::cout << "2-loop \\alpha^2      Mh= " << JKVH80.MH()  << ", [mW/MW -1] = " 
                << pow(alpha0/4./Pi,2)*W_ACOV_80.m20(0,0) << std::endl;
      std::cout << "2-loop \\alpha^2      Mh= " << JKVH200.MH()  << ", [mW/MW -1] = " 
                << pow(alpha0/4./Pi,2)*W_ACOV_200.m20(0,0) << std::endl;

      // 
      // Plot 2-loop bosonic part
      // 
      Plot1 plot("mW_0209084_fig6", "\\alpha^2 conversion between \\bar{MS} and  On-Shell W mass", "mH", "mmW(mZ)/MMW-1", "2-loop ");
      AlphaS as;
      long double mHstep  = 10; // GeV
      long double mHstart = 80; // GeV

      for (int mHi = 0; mHi < 13; mHi++)
        {
          // from Table.1
          SMinput DS2l(0, 80.419, 91.188, mHstart + mHi*mHstep, 174.3);
          
          WW dW  = WW(DS2l, DS2l.MMZ());          
          
          plot.add(DS2l.MH(), 
                   // alphaMt/4./Pi*dW.m10().real(), 
                   // alphaMt/4./Pi*as(DS2l.MMt())/4./Pi*dH.m11().real(),
                   pow(alpha0/4./Pi,2)*dW.m20(0,0).real());
        }


    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

