#include <iostream>
#include "mr.hpp"
#include "CRunDec.h"

#include "betaQEDQCD.hpp"


long double al(long double alMZ, long double MMZ, long double ass, long double mm, size_t Nu = 3, size_t Nd = 3, size_t Nl = 3)
{

  long double nu = static_cast<long double>(Nu);
  long double nd = static_cast<long double>(Nd);
  long double nl = static_cast<long double>(Nl);

  
  return alMZ/(1 - alMZ/4./Pi*(
                             7. - 4./3.*(3.*(nd + 4.*nu)/9. + nl)// -11./3.
                             -16.*(nd+4.*nu)/9.// -80./3.
                             *ass/4./Pi)*log(MMZ/mm));
}


int main (int argc, char *argv[])
{
  try
    {

      loglevel = logERROR;
      
      // Input: Pole masses and Fermi constant in OS scheme
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);
      
      // Running QCD coupling for as(Mt) from as(MZ)
      AlphaS as(oi);
      AlphaSolve  al(oi, pdg2014::Gf, as(oi.Mt()));      
      
      // Set of all running parameters at scale Mt
      P2MS pMSmt(oi,pdg2014::Gf, as(oi.Mt()), oi.Mt(), order::all);
      
      // Set of all running parameters at scale MZ
      P2MS pMSmZ(oi,pdg2014::Gf, as(oi.MZ()), oi.MZ(), order::all);

      // PDG2012 input for MH values from table ...
      std::vector<OSinput> oiMH;
      oiMH.push_back(oi.setMH(124));
      oiMH.push_back(oi.setMH(125));
      oiMH.push_back(oi.setMH(126));

      // And for physical value MH=125.66
      OSinput oiMHphys = oi.setMH(125.66);


      // AlphaS as;
      // CRunDec asRD;      
      
      // long double asMZ = pdg2014::asMZ;
      // long double aMZ  = pdg2014::aMZ;

      long double aMt  = al(oi.Mt());
      long double asMt = as(oi.Mt());
      
      
      // PDG2012

      std::cout << " Alpha at mu=Mt " << 1./aMt << std::endl;
      for (std::vector<OSinput>::iterator it = oiMH.begin(); it != oiMH.end(); ++it)
        {
          tt<OS> dMt(*it, it->MMt());
          
          P2MS pMSmtMH(*it,pdg2014::Gf, as(it->Mt()), it->Mt(), order::all);
          
          
          // Table
          std::cout << " | " << it->MH()
                    << " | " << it->Mt()*pow(asMt/4./Pi,1)*dMt.x01() + it->Mt()*pow(asMt/4./Pi,2)*dMt.x02() + it->Mt()*pow(asMt/4./Pi,3)*dMt.x03()
                    << " | " << it->Mt()*aMt/4./Pi*dMt.x10() 
                    << " | " << it->Mt()*aMt/4./Pi*asMt/4./Pi*dMt.x11() 
                    << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.xgl20() 
                    << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.x20() 
                    << std::endl;
          
        }

//       // PDG 2014
//       std::vector<OSinput> KPV;


//       KPV.push_back(OSinput(kmb*4.9, pdg2014::MW, pdg2014::MZ, 124, pdg2014::Mt));
//       KPV.push_back(OSinput(kmb*4.9, pdg2014::MW, pdg2014::MZ, 125, pdg2014::Mt));
//       KPV.push_back(OSinput(kmb*4.9, pdg2014::MW, pdg2014::MZ, 126, pdg2014::Mt));

//       OSinput KPVphys = OSinput(4.9, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);


//       // PDG2014
//       asRD.nfMmu[0].nf = 6;
//       asRD.nfMmu[0].Mth = KPVphys.Mt();
//       asRD.nfMmu[0].muth = KPVphys.Mt();
      
            
//       asMt = asRD.AlL2AlH(pdg2014::asMZ, KPVphys.MZ(), asRD.nfMmu, KPVphys.Mt(), 4);


//       std::pair<double,double> a2 = runQEDQCD(pdg2014::aMZ, pdg2014::asMZ, pdg2014::MZ, pdg2014::Mt);

//       AlphaQEDQCD FinAl(KPVphys);

//       std::cout << " as = " << a2.first << "     1/a = " << 1./a2.second << std::endl;
//       long double aMt2014 = al(pdg2014::aMZ, pow(pdg2014::MZ,2), pdg2014::asMZ, pow(pdg2014::Mt,2));
//       std::cout << "\\alpha              = " << 1./aMt2014 << std::endl;
//       std::cout << "\\alpha MZ           = " << 1./pdg2014::aMZ << std::endl;

//       std::cout << std::setprecision(10);
//       std::cout << "\\alpha MZ Final     = " << 1./FinAl.QED(pdg2014::Mt) << std::endl;
      
      
//       // AlphaQEDQCD aEMQCD(KPVphys);

      

//       // return 0;
      
//       AlphaS asPik(KPVphys); 

//       std::cout << "PDG2014" << std::endl;
//       std::cout << "\\aRD_S  (\\mu = Mt) = " << asMt  << std::endl;
//       std::cout << "\\aPIK_S  (\\mu = Mt) = " << asPik(pdg2014::Mt)  << std::endl;
      
      

//       // aMt = al(aMZ, KPVphys.MMZ(), pdg2014::asMZ, KPVphys.MMt());

//       aMt = 1./127.72063;       // GAPP value

//       aMt = FinAl.QED(pdg2014::Mt);
//       std::cout << "\\alpha              = " << 1./aMt << std::endl;
//       std::cout << "\\alpha00              = " << 1./pdg2014::aMZ << std::endl;
//       std::cout << "\\alpha10              = " << 1./al(pdg2014::aMZ, KPVphys.MMZ(), 0*pdg2014::asMZ, KPVphys.MMt()) << std::endl;
//       std::cout << "\\alpha11              = " << 1./al(pdg2014::aMZ, KPVphys.MMZ(), pdg2014::asMZ, KPVphys.MMt()) << std::endl;
      
//       for (std::vector<OSinput>::iterator it = KPV.begin(); it != KPV.end(); ++it)
//         {
//           tt<OS> dMt(*it, it->MMt());

//           dr<OS> dR(*it, it->MMt());

//           // Table
//           std::cout << " Mt-mt = | " << it->MH()
//                     << " | " << it->Mt()*pow(asMt/4./Pi,1)*dMt.x01() + it->Mt()*pow(asMt/4./Pi,2)*dMt.x02() + it->Mt()*pow(asMt/4./Pi,3)*dMt.x03()
//                     << " | " << it->Mt()*aMt/4./Pi*dMt.x10() 
//                     << " | " << it->Mt()*aMt/4./Pi*asMt/4./Pi*dMt.x11() 
//                     << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.xgl20() 
//                     << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.x20() 
//                     << " | " << it->Mt()*pow(aMt/4./Pi,2)*dMt.x20() 
// + it->Mt()*aMt/4./Pi*asMt/4./Pi*dMt.x11() 
// + it->Mt()*aMt/4./Pi*dMt.x10() 
// + it->Mt()*pow(asMt/4./Pi,1)*dMt.x01() + it->Mt()*pow(asMt/4./Pi,2)*dMt.x02() + it->Mt()*pow(asMt/4./Pi,3)*dMt.x03()
//                     << std::endl;

//           std::cout << " dr = | " << it->MH()
//                     << " | " << pow(aMt/4./Pi,2)*dR.drgl20() 
//                     << " | " << pow(aMt/4./Pi,2)*dR.dr20() 
//                     << std::endl;

//           std::cout << " Yt-yt *10^4= | " << it->MH()
//                     << " | " << 10000*(pow(asMt/4./Pi,1)*dMt.x01() + pow(asMt/4./Pi,2)*dMt.x02() + pow(asMt/4./Pi,3)*dMt.x03())
//                     << " | " << 10000*(aMt/4./Pi*dMt.y10() )
//                     << " | " << 10000*(aMt/4./Pi*asMt/4./Pi*dMt.y11()) 
//                     << " | " << 10000*(pow(aMt/4./Pi,2)*dMt.ygl20()) 
//                     << " | " << 10000*(pow(aMt/4./Pi,2)*dMt.y20() )
//                     << " | " << 10000*(pow(aMt/4./Pi,2)*dMt.y20()
// + aMt/4./Pi*asMt/4./Pi*dMt.y11() 
// + aMt/4./Pi*dMt.y10() 
// + pow(asMt/4./Pi,1)*dMt.x01() + pow(asMt/4./Pi,2)*dMt.x02() + pow(asMt/4./Pi,3)*dMt.x03())
//                     << std::endl;
          
//         }


//       std::cout << "\n\n Higgs: " << std::endl;
//       for (std::vector<OSinput>::iterator it = KPV.begin(); it != KPV.end(); ++it)
//         {
//           HH<OS> dMH  = HH<OS>(*it, it->MMt());
          
          
          
//           // Table
//           std::cout << " MH-mH = | " << it->MH()
//             // << " | " << it->Mt()*pow(asMt/4./Pi,1)*dMH.x01() + it->Mt()*pow(asMt/4./Pi,2)*dMH.x02() + it->Mt()*pow(asMt/4./Pi,3)*dMH.x03()
//                     << " | " << it->MH()*(sqrt(1 + aMt/4./Pi*dMH.x10()) - 1)
//                     << " | " << it->MH()*(sqrt(1 + aMt/4./Pi*asMt/4./Pi*dMH.x11()) - 1) 
//                     << " | " << it->MH()*(sqrt(1 + pow(aMt/4./Pi,2)*dMH.xgl20()) - 1) 
//                     << " | " << it->MH()*(sqrt(1 + pow(aMt/4./Pi,2)*dMH.x20()) - 1) 
//                     << " | " << it->MH()*(sqrt(1 + pow(aMt/4./Pi,2)*dMH.x20() 
//                                                + aMt/4./Pi*asMt/4./Pi*dMH.x11() 
//                                                + aMt/4./Pi*dMH.x10() ) - 1)
//             // + it->Mt()*pow(asMt/4./Pi,1)*dMH.x01() + it->Mt()*pow(asMt/4./Pi,2)*dMH.x02() + it->Mt()*pow(asMt/4./Pi,3)*dMH.x03()
//                     << std::endl;
          
//           // Table
//           std::cout << " MH^2-mH^2 = | " << it->MH()
//             // << " | " << it->Mt()*pow(asMt/4./Pi,1)*dMH.x01() + it->Mt()*pow(asMt/4./Pi,2)*dMH.x02() + it->Mt()*pow(asMt/4./Pi,3)*dMH.x03()
//                     << " | " << it->MMH()*aMt/4./Pi*dMH.x10()
//                     << " | " << it->MMH()*aMt/4./Pi*asMt/4./Pi*dMH.x11()
//                     << " | " << it->MMH()*pow(aMt/4./Pi,2)*dMH.xgl20()
//                     << " | " << it->MMH()*pow(aMt/4./Pi,2)*dMH.x20() 
//                     << " | " << it->MMH()*(pow(aMt/4./Pi,2)*dMH.x20() 
//                                            + aMt/4./Pi*asMt/4./Pi*dMH.x11() 
//                                            + aMt/4./Pi*dMH.x10())
//             // + it->Mt()*pow(asMt/4./Pi,1)*dMH.x01() + it->Mt()*pow(asMt/4./Pi,2)*dMH.x02() + it->Mt()*pow(asMt/4./Pi,3)*dMH.x03()
//                     << std::endl;
          
          
//           std::cout << " Lam-lam = | " << it->MH()
//             // << " | " << pow(asMt/4./Pi,1)*dMH.x01() + pow(asMt/4./Pi,2)*dMH.x02() + pow(asMt/4./Pi,3)*dMH.x03()
//                     << " | " << 10000*(aMt/4./Pi*dMH.y10()) 
//                     << " | " << 10000*(aMt/4./Pi*asMt/4./Pi*dMH.y11() )
//                     << " | " << 10000*(pow(aMt/4./Pi,2)*dMH.ygl20() )
//                     << " | " << 10000*(pow(aMt/4./Pi,2)*dMH.y20() )
//                     << " | " << 10000*(pow(aMt/4./Pi,2)*dMH.y20() 
//                                        + aMt/4./Pi*asMt/4./Pi*dMH.y11() 
//                                        + aMt/4./Pi*dMH.y10()) 
//             // + pow(asMt/4./Pi,1)*dMH.x01() + pow(asMt/4./Pi,2)*dMH.x02() + pow(asMt/4./Pi,3)*dMH.x03()
//                     << std::endl;
          
//         }
      
      
//       bb<OS> dMb(KPVphys, KPVphys.MMb());
      
      
//       // PDG2014
//       // asRD.nfMmu[0].nf = 5;
//       // asRD.nfMmu[0].Mth = KPVphys.Mt();
//       // asRD.nfMmu[0].muth = KPVphys.Mt();
      
//       CRunDec asRD5(5);      
      
      
//       long double asMb = asRD5.AlphasExact(pdg2014::asMZ, KPVphys.MZ(), KPVphys.Mb(), 4);

//       AlphaS  asMBPIK(KPVphys);

//       std::cout << "PDG2014" << std::endl;
//       std::cout << "\\aRD_S  (\\mu = Mb) = " << asMb  << std::endl;
//       std::cout << "\\aPIK_S  (\\mu = Mb) = " << asMBPIK(4.9)  << std::endl;
      
//       long double aMb =       aMt = FinAl.QED(pdg2014::Mb);
//       // al(aMZ, KPVphys.MMZ(), asMb, KPVphys.MMb());
//       std::cout << "\\alpha a(MB)             = " << 1./aMb << std::endl;
      
//       long double aMbMZ = al(aMZ, KPVphys.MMZ(), pdg2014::asMZ, KPVphys.MMb());
//       std::cout << "\\alpha as(MZ)             = " << 1./aMbMZ << std::endl;

      
//       // Table
//       std::cout << " Mb-mb = | " << KPVphys.Mb()
//                 << " | " << KPVphys.Mb()*pow(asMb/4./Pi,1)*dMb.x01() + KPVphys.Mb()*pow(asMb/4./Pi,2)*dMb.x02() + KPVphys.Mb()*pow(asMb/4./Pi,3)*dMb.x03()
//                 << " | " << KPVphys.Mb()*aMb/4./Pi*dMb.x10() 
//                 << " | " << KPVphys.Mb()*aMb/4./Pi*asMb/4./Pi*dMb.x11() 
//                 << " | " << KPVphys.Mb()*pow(aMb/4./Pi,2)*dMb.xgl20() 
//                 << " | " << KPVphys.Mb()*pow(aMb/4./Pi,2)*dMb.x20() 
//                 << " | " << KPVphys.Mb()*pow(aMb/4./Pi,2)*dMb.x20() 
//         + KPVphys.Mb()*aMb/4./Pi*asMb/4./Pi*dMb.x11() 
//         + KPVphys.Mb()*aMb/4./Pi*dMb.x10() 
//         + KPVphys.Mb()*pow(asMb/4./Pi,1)*dMb.x01() + KPVphys.Mb()*pow(asMb/4./Pi,2)*dMb.x02() + KPVphys.Mb()*pow(asMb/4./Pi,3)*dMb.x03()
//                     << std::endl;


//       std::cout << " Yb-yb = | "
//                 << " | " << pow(asMb/4./Pi,1)*dMb.y01() + pow(asMb/4./Pi,2)*dMb.x02() + pow(asMb/4./Pi,3)*dMb.x03()
//                 << " | " << aMb/4./Pi*dMb.y10() 
//                 << " | " << aMb/4./Pi*asMb/4./Pi*dMb.y11() 
//                 << " | " << pow(aMb/4./Pi,2)*dMb.ygl20() 
//                 << " | " << pow(aMb/4./Pi,2)*dMb.y20() 
//                 << " | " << pow(aMb/4./Pi,2)*dMb.y20() 
//         + aMb/4./Pi*asMb/4./Pi*dMb.y11() 
//         + aMb/4./Pi*dMb.y10() 
//         + pow(asMb/4./Pi,1)*dMb.y01() + pow(asMb/4./Pi,2)*dMb.x02() + pow(asMb/4./Pi,3)*dMb.x03()
//                     << std::endl;


      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

