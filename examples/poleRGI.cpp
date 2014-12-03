#include "mr.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);

      // Input from Martin
      long double Mt = 173.10;
      MSinput SPM = MSinput::fromConsts( pow(Mt,2), // Input scale
                                         sqrt(2.)*93.36,
                                         0.127,
                                         0, 0.936, 0.648, 0.358);

      // 0.127 247.0 0.936 1.167 0.648 0.358 173.1

      std::cout << "Mb= " << SPM.mb() << std::endl; 
      std::cout << "MW= " << SPM.mW() << std::endl; 
      std::cout << "MZ= " << SPM.mZ() << std::endl; 
      std::cout << "MH= " << SPM.mH() << std::endl; 
      std::cout << "Mt= " << SPM.mt() << std::endl;
      std::cout << "V=  " << SPM.vev() << std::endl; 
      std::cout << "1/al= " << 1./SPM.alpha() << std::endl; 


      AlphaS as(pdg2014::MZ, pdg2014::asMZ, 4, 0, Mt);

      HH<MS> H_mt(SPM, pow(Mt,2));
      
      
      double MH10,MH11,MH20;

      MH10 = sqrt(SPM.mmH() * (1 +
                               SPM.alpha()/4./Pi*H_mt.x10()
                               ).real());
      
      MH11 = sqrt(SPM.mmH() * (1 +
                               SPM.alpha()/4./Pi*H_mt.x10() +
                               SPM.alpha()/4./Pi*as(Mt)/4./Pi*H_mt.x11()
                               ).real());
      
      MH20 = sqrt(SPM.mmH() * (1 +
                               SPM.alpha()/4./Pi*H_mt.x10() +
                               SPM.alpha()/4./Pi*as(Mt)/4./Pi*H_mt.x11() +
                               pow(SPM.alpha()/4./Pi,2)*H_mt.x20()
                               ).real());
      
      std::cout << "\n MH[ EH ]     = " << MH10 << std::endl;
      std::cout << " MW[ EW*QCD ] = " << MH11 << std::endl;
      std::cout << " MH[ EW*EW ]  = " << MH20 << std::endl;
      
      
      return 0;

      
      double muIn = 91.1876;

      
      CouplingsMu av(
                     5./3.*pow(0.3497/4./Pi,2), // GUT normalization
                     pow(0.6530/4./Pi,2),
                     pow(1.2200/4./Pi,2),
                     pow(0.9347/4./Pi,2),
                     pow(0.0238/4./Pi,2),
                     pow(0.0104/4./Pi,2),
                     0.8070/6.*pow(4.*Pi,-2),
                     // no vev in input
                     89.096*sqrt(2.),
                     pow(muIn,2),
                     3);
      

      


      state_type avZ = av(pow(muIn,2));

      std::cout << 4*Pi*sqrt(3./5.*avZ[0]) << std::endl;
      std::cout << avZ[7] << std::endl;
      
      MSinput msi = MSinput::fromConsts( pow(muIn,2), // Input scale
                                         avZ[7],
                                         pow(4.*Pi,2)*avZ[6],
                                         0,
                                         4*Pi*sqrt(avZ[3]),
                                         4*Pi*sqrt(avZ[1]),
                                         4*Pi*sqrt(3./5.*avZ[0]));

      std::cout << "Mb= " << msi.mb() << std::endl; 
      std::cout << "MW= " << msi.mW() << std::endl; 
      std::cout << "MZ= " << msi.mZ() << std::endl; 
      std::cout << "MH= " << msi.mH() << std::endl; 
      std::cout << "Mt= " << msi.mt() << std::endl; 
      std::cout << "1/al= " << 1./msi.alpha() << std::endl; 


      // AlphaS as(pdg2014::MZ, pdg2014::asMZ, 4, 0, pdg2014::Mt);

      WW<MS> W_mt(msi, pow(muIn,2));
      
      std::cout << " MW = " << sqrt(msi.mmW() * (1 +
                                                 msi.alpha()/4./Pi*W_mt.x10() +
                                                 msi.alpha()/4./Pi*as(muIn)/4./Pi*W_mt.x11() +
                                                 pow(msi.alpha()/4./Pi,2)*W_mt.x20()
                                                 ).real())
                << std::endl;


      //
      //
      //
      //
      //
      // 

      Plot3 plotMW("rgeMW", "\\alpha*\\alpha_S conversion between \\bar{MS} and  On-Shell Higgs mass", "Q", "MH", "1-loop ", "2-loop EW*QCD", "2-loop EW^2");
     
      for (int mu = 100; mu <= 250; mu += 10)
        {
          
          double muOut = double(mu);
          
          state_type avMU = av(pow(muOut,2));
          
          std::cout << 4*Pi*sqrt(3./5.*avMU[0]) << std::endl;
          std::cout << avMU[7] << std::endl;
          
          MSinput msMU = MSinput::fromConsts( pow(muOut,2), // Input scale
                                              avMU[7],
                                              pow(4.*Pi,2)*avMU[6],
                                              0,
                                              4*Pi*sqrt(avMU[3]),
                                              4*Pi*sqrt(avMU[1]),
                                              4*Pi*sqrt(3./5.*avMU[0]));
          
          // std::cout << "Mb= " << msMU.mb() << std::endl; 
          // std::cout << "MW= " << msMU.mW() << std::endl; 
          // std::cout << "MZ= " << msMU.mZ() << std::endl; 
          // std::cout << "MH= " << msMU.mH() << std::endl; 
          // std::cout << "Mt= " << msMU.mt() << std::endl; 
          // std::cout << "1/al= " << 1./msMU.alpha() << std::endl; 
          
          
          // AlphaS as(pdg2014::MZ, pdg2014::asMZ, 4, 0, pdg2014::Mt);
          
          WW<MS> W_mMU(msMU, pow(muOut,2));
          double MW10,MW11,MW20;

          MW10 = sqrt(msMU.mmW() * (1 +
                                    msMU.alpha()/4./Pi*W_mMU.x10()
                                    ).real());
          
          MW11 = sqrt(msMU.mmW() * (1 +
                                    msMU.alpha()/4./Pi*W_mMU.x10() +
                                    msMU.alpha()/4./Pi*as(muOut)/4./Pi*W_mMU.x11()
                                    ).real());

          MW20 = sqrt(msMU.mmW() * (1 +
                                    msMU.alpha()/4./Pi*W_mMU.x10() +
                                    msMU.alpha()/4./Pi*as(muOut)/4./Pi*W_mMU.x11() +
                                    pow(msMU.alpha()/4./Pi,2)*W_mMU.x20()
                                    ).real());
          
          std::cout << "\n MW[ EW ]     = " << MW10 << std::endl;
          std::cout << " MW[ EW*QCD ] = " << MW11 << std::endl;
          std::cout << " MW[ EW*EW ]  = " << MW20 << std::endl;

          plotMW.add(muOut,MW10,MW11,MW20);
        } 
        }
      catch (std::exception &p) 
        {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

