#include <iostream>
#include "mr.hpp"
#include "CRunDec.h"

/*#define asMz 0.1184
#define Mz   91.1876
#define Mt   173.5
#define Mb   4.4
*/
// #define Mc   1.5
// #define muc  1.279
// #define mub  4.163
// #define Mtau 1.777


long double al(long double alMZ, long double MMZ, long double ass, long double mm)
{
  return alMZ/(1-alMZ/4./Pi*(-11./3.-80./3.*ass/4./Pi)*log(MMZ/mm));
}


int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);

      long double kmb = 1.;

      // Compare with:
      std::vector<OSinput> KV;
      KV.push_back(OSinput(kmb*4.4, 80.385, 91.1876, 124, 173.5));
      KV.push_back(OSinput(kmb*4.4, 80.385, 91.1876, 125, 173.5));
      KV.push_back(OSinput(kmb*4.4, 80.385, 91.1876, 126, 173.5));

      OSinput KVphys = OSinput(4.4, 80.385, 91.1876, 125.6, 173.5);


      AlphaS as;
      CRunDec asRD;      
      
      long double asMZ = 0.1184;
      long double aMZ  = 1./127.944;

      long double aMt = 1./127.944;
      long double asMt;
      
      
      // PDG2012
      asRD.nfMmu[0].nf = 6;
      asRD.nfMmu[0].Mth = KVphys.Mt();
      asRD.nfMmu[0].muth = KVphys.Mt();
      
            
      asMt = asRD.AlL2AlH(asMZ, KVphys.MZ(), asRD.nfMmu, KVphys.Mt(), 4);
      std::cout << "PDG2012" << std::endl;
      std::cout << "\\aRD_S  (\\mu = Mt) = " << asMt  << std::endl;
      
      aMt = al(aMZ, KVphys.MMZ(), asMZ, KVphys.MMt());
      std::cout << "\\alpha              = " << 1./aMt << std::endl;
      
      for (std::vector<OSinput>::iterator it = KV.begin(); it != KV.end(); ++it)
        {
          tt dMt  = tt(*it, it->MMt());



          // Table
          std::cout << " | " << it->MH()
                    << " | " << it->Mt()*pow(asMt/4./Pi,1)*dMt.m01().real() + it->Mt()*pow(asMt/4./Pi,2)*dMt.m02().real() + it->Mt()*pow(asMt/4./Pi,3)*dMt.m03().real()
                    << " | " << it->Mt()*aMt/4./Pi*dMt.m10().real() 
                    << " | " << it->Mt()*aMt/4./Pi*asMt/4./Pi*dMt.m11().real() 
                    << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.mgl20().real() 
                    << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.m20().real() 
                    << std::endl;
          
        }

      // PDG 2014
      std::vector<OSinput> KPV;


      KPV.push_back(OSinput(kmb*4.9, 80.385, 91.1876, 124, 173.2));
      KPV.push_back(OSinput(kmb*4.9, 80.385, 91.1876, 125, 173.2));
      KPV.push_back(OSinput(kmb*4.9, 80.385, 91.1876, 126, 173.2));

      OSinput KPVphys = OSinput(4.9, 80.385, 91.1876, 125.7, 173.2);


      // PDG2014
      asRD.nfMmu[0].nf = 6;
      asRD.nfMmu[0].Mth = KPVphys.Mt();
      asRD.nfMmu[0].muth = KPVphys.Mt();
      
            
      asMt = asRD.AlL2AlH(asMZ, KPVphys.MZ(), asRD.nfMmu, KPVphys.Mt(), 4);

      std::cout << "PDG2014" << std::endl;
      std::cout << "\\aRD_S  (\\mu = Mt) = " << asMt  << std::endl;
      

      aMt = al(aMZ, KPVphys.MMZ(), asMZ, KPVphys.MMt());

      std::cout << "\\alpha              = " << 1./aMt << std::endl;

      for (std::vector<OSinput>::iterator it = KPV.begin(); it != KPV.end(); ++it)
        {
          tt dMt  = tt(*it, it->MMt());



          // Table
          std::cout << " Mt-mt = | " << it->MH()
                    << " | " << it->Mt()*pow(asMt/4./Pi,1)*dMt.m01().real() + it->Mt()*pow(asMt/4./Pi,2)*dMt.m02().real() + it->Mt()*pow(asMt/4./Pi,3)*dMt.m03().real()
                    << " | " << it->Mt()*aMt/4./Pi*dMt.m10().real() 
                    << " | " << it->Mt()*aMt/4./Pi*asMt/4./Pi*dMt.m11().real() 
                    << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.mgl20().real() 
                    << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.m20().real() 
                    << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.m20().real() 
+ it->Mt()*aMt/4./Pi*asMt/4./Pi*dMt.m11().real() 
+ it->Mt()*aMt/4./Pi*dMt.m10().real() 
+ it->Mt()*pow(asMt/4./Pi,1)*dMt.m01().real() + it->Mt()*pow(asMt/4./Pi,2)*dMt.m02().real() + it->Mt()*pow(asMt/4./Pi,3)*dMt.m03().real()
                    << std::endl;


          std::cout << " Yt-yt *10^4= | " << it->MH()
                    << " | " << 10000*(pow(asMt/4./Pi,1)*dMt.m01().real() + pow(asMt/4./Pi,2)*dMt.m02().real() + pow(asMt/4./Pi,3)*dMt.m03().real())
                    << " | " << 10000*(aMt/4./Pi*dMt.my10().real() )
                    << " | " << 10000*(aMt/4./Pi*asMt/4./Pi*dMt.my11().real()) 
                    << " | " << 10000*(pow(aMt/4./Pi,2)*dMt.mygl20().real()) 
                    << " | " << 10000*(pow(aMt/4./Pi,2)*dMt.my20().real() )
                    << " | " << 10000*(pow(aMt/4./Pi,2)*dMt.my20().real()
+ aMt/4./Pi*asMt/4./Pi*dMt.my11().real() 
+ aMt/4./Pi*dMt.my10().real() 
+ pow(asMt/4./Pi,1)*dMt.m01().real() + pow(asMt/4./Pi,2)*dMt.m02().real() + pow(asMt/4./Pi,3)*dMt.m03().real())
                    << std::endl;
          
        }


std::cout << "\n\n Higgs: " << std::endl;
      for (std::vector<OSinput>::iterator it = KPV.begin(); it != KPV.end(); ++it)
        {
          HH<OS> dMH  = HH<OS>(*it, it->MMt());



          // Table
          std::cout << " MH-mH = | " << it->MH()
                    // << " | " << it->Mt()*pow(asMt/4./Pi,1)*dMH.m01().real() + it->Mt()*pow(asMt/4./Pi,2)*dMH.m02().real() + it->Mt()*pow(asMt/4./Pi,3)*dMH.m03().real()
                    << " | " << it->MH()*(sqrt(1 + aMt/4./Pi*dMH.m10().real()) - 1)
                    << " | " << it->MH()*(sqrt(1 + aMt/4./Pi*asMt/4./Pi*dMH.m11().real()) - 1) 
                    << " | " << it->MH()*(sqrt(1 + pow(aMt/4./Pi,2)*dMH.mgl20().real()) - 1) 
                    << " | " << it->MH()*(sqrt(1 + pow(aMt/4./Pi,2)*dMH.m20().real()) - 1) 
                    << " | " << it->MH()*(sqrt(1 + pow(aMt/4./Pi,2)*dMH.m20().real() 
+ aMt/4./Pi*asMt/4./Pi*dMH.m11().real() 
+ aMt/4./Pi*dMH.m10().real() ) - 1)
// + it->Mt()*pow(asMt/4./Pi,1)*dMH.m01().real() + it->Mt()*pow(asMt/4./Pi,2)*dMH.m02().real() + it->Mt()*pow(asMt/4./Pi,3)*dMH.m03().real()
                    << std::endl;

          // Table
          std::cout << " MH^2-mH^2 = | " << it->MH()
                    // << " | " << it->Mt()*pow(asMt/4./Pi,1)*dMH.m01().real() + it->Mt()*pow(asMt/4./Pi,2)*dMH.m02().real() + it->Mt()*pow(asMt/4./Pi,3)*dMH.m03().real()
                    << " | " << it->MMH()*aMt/4./Pi*dMH.m10().real()
                    << " | " << it->MMH()*aMt/4./Pi*asMt/4./Pi*dMH.m11().real()
                    << " | " << it->MMH()*pow(aMt/4./Pi,2)*dMH.mgl20().real()
                    << " | " << it->MMH()*pow(aMt/4./Pi,2)*dMH.m20().real() 
                    << " | " << it->MMH()*(pow(aMt/4./Pi,2)*dMH.m20().real() 
                                           + aMt/4./Pi*asMt/4./Pi*dMH.m11().real() 
                                           + aMt/4./Pi*dMH.m10().real())
// + it->Mt()*pow(asMt/4./Pi,1)*dMH.m01().real() + it->Mt()*pow(asMt/4./Pi,2)*dMH.m02().real() + it->Mt()*pow(asMt/4./Pi,3)*dMH.m03().real()
                    << std::endl;


          std::cout << " Lam-lam = | " << it->MH()
                    // << " | " << pow(asMt/4./Pi,1)*dMH.m01().real() + pow(asMt/4./Pi,2)*dMH.m02().real() + pow(asMt/4./Pi,3)*dMH.m03().real()
                    << " | " << 10000*(aMt/4./Pi*dMH.my10().real()) 
                    << " | " << 10000*(aMt/4./Pi*asMt/4./Pi*dMH.my11().real() )
                    << " | " << 10000*(pow(aMt/4./Pi,2)*dMH.mygl20().real() )
                    << " | " << 10000*(pow(aMt/4./Pi,2)*dMH.my20().real() )
                    << " | " << 10000*(pow(aMt/4./Pi,2)*dMH.my20().real() 
                                       + aMt/4./Pi*asMt/4./Pi*dMH.my11().real() 
                                       + aMt/4./Pi*dMH.my10().real()) 
            // + pow(asMt/4./Pi,1)*dMH.m01().real() + pow(asMt/4./Pi,2)*dMH.m02().real() + pow(asMt/4./Pi,3)*dMH.m03().real()
                    << std::endl;
          
        }


      bb dMb  = bb(KPVphys, KPVphys.MMb());
      
      
      // PDG2014
      // asRD.nfMmu[0].nf = 5;
      // asRD.nfMmu[0].Mth = KPVphys.Mt();
      // asRD.nfMmu[0].muth = KPVphys.Mt();
      
      CRunDec asRD5(5);      
      
      
            
      long double asMb = asRD5.AlphasExact(asMZ, KPVphys.MZ(), KPVphys.Mb(), 4);

      std::cout << "PDG2014" << std::endl;
      std::cout << "\\aRD_S  (\\mu = Mb) = " << asMb  << std::endl;
      
      long double aMb = al(aMZ, KPVphys.MMZ(), asMb, KPVphys.MMb());
      std::cout << "\\alpha              = " << 1./aMb << std::endl;

      
      // Table
      std::cout << " Mb-mb = | " << KPVphys.Mb()
                << " | " << KPVphys.Mb()*pow(asMb/4./Pi,1)*dMb.m01().real() + KPVphys.Mb()*pow(asMb/4./Pi,2)*dMb.m02().real() + KPVphys.Mb()*pow(asMb/4./Pi,3)*dMb.m03().real()
                << " | " << KPVphys.Mb()*aMb/4./Pi*dMb.m10().real() 
                << " | " << KPVphys.Mb()*aMb/4./Pi*asMb/4./Pi*dMb.m11().real() 
                << " | " << KPVphys.Mb()*pow(aMb/4./Pi,2)*dMb.mgl20().real() 
                << " | " << KPVphys.Mb()*pow(aMb/4./Pi,2)*dMb.m20().real() 
                << " | " << KPVphys.Mb()*pow(aMb/4./Pi,2)*dMb.m20().real() 
        + KPVphys.Mb()*aMb/4./Pi*asMb/4./Pi*dMb.m11().real() 
        + KPVphys.Mb()*aMb/4./Pi*dMb.m10().real() 
        + KPVphys.Mb()*pow(asMb/4./Pi,1)*dMb.m01().real() + KPVphys.Mb()*pow(asMb/4./Pi,2)*dMb.m02().real() + KPVphys.Mb()*pow(asMb/4./Pi,3)*dMb.m03().real()
                    << std::endl;


      std::cout << " Yb-yb = | "
                << " | " << pow(asMb/4./Pi,1)*dMb.my01().real() + pow(asMb/4./Pi,2)*dMb.m02().real() + pow(asMb/4./Pi,3)*dMb.m03().real()
                << " | " << aMb/4./Pi*dMb.my10().real() 
                << " | " << aMb/4./Pi*asMb/4./Pi*dMb.my11().real() 
                << " | " << pow(aMb/4./Pi,2)*dMb.mygl20().real() 
                << " | " << pow(aMb/4./Pi,2)*dMb.my20().real() 
                << " | " << pow(aMb/4./Pi,2)*dMb.my20().real() 
        + aMb/4./Pi*asMb/4./Pi*dMb.my11().real() 
        + aMb/4./Pi*dMb.my10().real() 
        + pow(asMb/4./Pi,1)*dMb.my01().real() + pow(asMb/4./Pi,2)*dMb.m02().real() + pow(asMb/4./Pi,3)*dMb.m03().real()
                    << std::endl;


      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

