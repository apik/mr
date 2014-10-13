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
      SMinput JKVH80(0, 80.419, 91.188, 80, 174.3);
      SMinput JKVH200(0, 80.419, 91.188, 200, 174.3);
         
      alphaMt  = 0.00779305;
      alphaS = 0.1184;

      long double alpha0 = 1./137.035999;

      ZZ Z80  = ZZ(JKVH80, JKVH80.MMZ());
      ZZ Z200 = ZZ(JKVH200, JKVH200.MMZ());

      std::cout << "\\mu = MW" << std::endl;
      // 
      // Bosonic part 1-loop
      // 
      std::cout << "No fermions: 1-loop \\alpha      Mh= " << JKVH80.MH()  << ", [mZ/MZ -1] = " 
                << alpha0/4./Pi*Z80.m10(0,0) << std::endl;
      std::cout << "No fermions: 1-loop \\alpha      Mh= " << JKVH80.MH()  << ", [mZ/MZ -1] = " 
                << alpha0/4./Pi*Z200.m10(0,0) << std::endl;

      // 
      // Full SM
      // 
      std::cout << "SM 1-loop \\alpha      Mh= " << JKVH80.MH()  << ", [mZ/MZ -1] = " 
                << alpha0/4./Pi*Z80.m10(2,1) << std::endl;
      std::cout << "SM 1-loop \\alpha      Mh= " << JKVH80.MH()  << ", [mZ/MZ -1] = " 
                << alpha0/4./Pi*Z200.m10(2,1) << std::endl;

      std::cout << "SM 2-loop \\alpha^2      Mh= " << JKVH80.MH()  << ", [mZ/MZ -1] = " 
                << pow(alpha0/4./Pi,2)*Z80.m20(0,0) << std::endl;
      std::cout << "SM 2-loop \\alpha^2      Mh= " << JKVH80.MH()  << ", [mZ/MZ -1] = " 
                << pow(alpha0/4./Pi,2)*Z200.m20(0,0) << std::endl;


      // In ACOV \mu^2 = MMZ
      ZZ Z_ACOV_80  = ZZ(JKVH80, JKVH80.MMZ());
      ZZ Z_ACOV_200 = ZZ(JKVH200, JKVH200.MMZ());
      
      // 
      // Bosonic part 2-loop, compare with ACOV: arXiv:hep-ph/0209084,
      // FIG. 6
      // 
      std::cout << "\\mu = MZ" << std::endl;
      std::cout << "2-loop \\alpha^2      Mh= " << JKVH80.MH()  << ", [mZ/MZ -1] = " 
                << pow(alpha0/4./Pi,2)*Z_ACOV_80.m20(0,0) << std::endl;
      std::cout << "2-loop \\alpha^2      Mh= " << JKVH200.MH()  << ", [mZ/MZ -1] = " 
                << pow(alpha0/4./Pi,2)*Z_ACOV_200.m20(0,0) << std::endl;

      // 
      // Plot 2-loop bosonic part
      // 
      Plot1 plot("mZ_0209084_fig6", "\\alpha^2 conversion between \\bar{MS} and  On-Shell Z mass", "mH", "mmZ(mZ)/MMZ-1", "2-loop ");
      AlphaS as;
      long double mHstep  = 10; // GeV
      long double mHstart = 80; // GeV

      for (int mHi = 0; mHi < 13; mHi++)
        {
          // from Table.1
          SMinput DS2l(0, 80.419, 91.188, mHstart + mHi*mHstep, 174.3);
          
          ZZ dZ  = ZZ(DS2l, DS2l.MMZ());          
          
          plot.add(DS2l.MH(), 
                   // alphaMt/4./Pi*dW.m10().real(), 
                   // alphaMt/4./Pi*as(DS2l.MMt())/4./Pi*dH.m11().real(),
                   pow(alpha0/4./Pi,2)*dZ.m20(0,0).real());
        }


      // 
      // Plot 2-loop bosonic part
      // 
      Plot4 plotJKV("mZ_0105304_fig3", "\\alpha^2 conversion between \\bar{MS} and  On-Shell Z mass", "mH", "mmZ(mZ)/MMZ-1", "1-loop nf=0", "1-loop nf=3", "2-loop", "1+2-loop");
      mHstep  = 100; // GeV
      mHstart = 100; // GeV
      
      for (int mHi = 0; mHi < 10; mHi++)
        {
          // from Table.1
          SMinput DS2l(4.4, 80.419, 91.188, mHstart + mHi*mHstep, 174.3);
          ZZ dZ  = ZZ(DS2l, DS2l.MMZ());          
          
          long double mZ_2l = pow(alpha0/4./Pi,2)*dZ.m20(0,0).real();
          long double mZ_1l = pow(alpha0/4./Pi,1)*dZ.m10(0,0).real();
          long double mZ_1lnf3 = pow(alpha0/4./Pi,1)*dZ.m10().real();

          std::cout << "mH = " << mHstart + mHi*mHstep << " mZ/MZ(MZ,nf=0)[1-loop] = " << mZ_1l << std::endl;
          std::cout << "mH = " << mHstart + mHi*mHstep << " mZ/MZ(MZ,nf=3)[1-loop] = " << mZ_1lnf3 << std::endl;
          std::cout << "mH = " << mHstart + mHi*mHstep << " mZ/MZ(MZ,nf=0)[2-loop] = " << mZ_2l << std::endl << std::endl;

          plotJKV.add(DS2l.MH(), 
                   mZ_1l,
                   mZ_1lnf3,
                   mZ_2l,
                   mZ_1l+mZ_2l
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

