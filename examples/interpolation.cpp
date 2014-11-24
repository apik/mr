// #include <CRunDec.h>
#include <iostream>
#include "mr.hpp"

#include <map>
#include <utility>

// Fit staff
#include <TGraph.h>
#include <TApplication.h>
#include <TFitResult.h>


int main (int argc, char *argv[])
{
  try
    {

      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS,alphaSMt;

      // Scale inv test
      long double a = 40.;
      // Compare with:
      alphaMt    = 1./127.72063;
      alphaSMt   = 0.10798036651966895;
      
      // Initial input
      // OSinput KPV = OSinput(4.4, 80.385, 91.1876, 125.6, 173.5);

      // PDG 2014
      OSinput KPV = OSinput(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);

      long double unitStep = 0.0001;
      long double dMtStep  = 0.51;
      long double dMHStep  = 0.4;
      long double dMWStep  = 0.015;


      // Every delta is 1 GeV
      // OSinput divMt = OSinput(4.4, 80.385, 91.1876, 125.6, 173.5 + dMtStep);
      // OSinput divMH = OSinput(4.4, 80.385, 91.1876, 125.6 + dMHStep, 173.5);
      // OSinput divMW = OSinput(4.4, 80.385 + dMWStep, 91.1876, 125.6, 173.5);

      OSinput divMt = OSinput(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt + dMtStep);
      OSinput divMH = OSinput(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH + dMHStep, pdg2014::Mt);
      OSinput divMW = OSinput(pdg2014::Mb, pdg2014::MW + dMWStep, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);

      PoleMass* xMbase0      = new tt(KPV  , KPV.MMZ());
      // tt xMttdivMt  = tt(divMt, KPV.MMZ());
      // tt xMttdivMH  = tt(divMH, KPV.MMZ());
      // tt xMttdivMW  = tt(divMW, KPV.MMZ());


      // Loop over energy scale points:

      long double mm[2] = {KPV.MMZ()};// , KPV.MMt()};
      
      typedef std::map<std::pair<std::string,long double>, std::vector<PoleMass*> > Mstore;
      Mstore mMap;
 
      for (int i = 0; i < 1; i++)
        {
           
          std::vector<PoleMass*> pmvWW(4);
          pmvWW[0] = new WW<OS>(KPV, mm[i]);
          pmvWW[1] = new WW<OS>(divMt, mm[i]);
          pmvWW[2] = new WW<OS>(divMH, mm[i]);
          pmvWW[3] = new WW<OS>(divMW, mm[i]);
          
          mMap[std::make_pair("W",mm[i])] = pmvWW;

          std::vector<PoleMass*> pmvZZ(4);
          pmvZZ[0] = new ZZ<OS>(KPV, mm[i]);
          pmvZZ[1] = new ZZ<OS>(divMt, mm[i]);
          pmvZZ[2] = new ZZ<OS>(divMH, mm[i]);
          pmvZZ[3] = new ZZ<OS>(divMW, mm[i]);
          
          mMap[std::make_pair("Z",mm[i])] = pmvZZ;

          std::vector<PoleMass*> pmvHH(4);
          pmvHH[0] = new HH<OS>(KPV, mm[i]);
          pmvHH[1] = new HH<OS>(divMt, mm[i]);
          pmvHH[2] = new HH<OS>(divMH, mm[i]);
          pmvHH[3] = new HH<OS>(divMW, mm[i]);
          
          mMap[std::make_pair("H",mm[i])] = pmvHH;

          std::vector<PoleMass*> pmvtt(4);
          pmvtt[0] = new tt(KPV, mm[i]);
          pmvtt[1] = new tt(divMt, mm[i]);
          pmvtt[2] = new tt(divMH, mm[i]);
          pmvtt[3] = new tt(divMW, mm[i]);
          
          mMap[std::make_pair("t",mm[i])] = pmvtt;
          
        }

      std::cout << " Terms = " << mMap.size() << std::endl;
      

      for (Mstore::const_iterator it = mMap.begin(); it != mMap.end(); ++it)
        for(size_t aspow = 0; aspow <= 1; aspow++)
          for(size_t apow = 2 - aspow; apow <= 2 - aspow ; apow++)
            {
              long double M0  = it->second[0]->x(apow,aspow);
              long double Mdt = it->second[1]->x(apow,aspow);
              long double MdH = it->second[2]->x(apow,aspow);
              long double MdW = it->second[3]->x(apow,aspow);
              // std::cout << "Interpolation for X" << it->first.first 
              //           << "["<< apow << "," << aspow <<"] at  scale \\mu = " 
              //           << sqrt(it->first.second) << std::endl;
              // std::cout << (Mdt - M0)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
              // std::cout << (MdH - M0)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
              // std::cout << (MdW - M0)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
              // std::cout << M0 << std::endl;      

              std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
              std::cout << "X_" << it->first.first 
                        << "^{"<< apow << "," << aspow <<"}(\\mu = " 
                        << sqrt(it->first.second) << ") & = & " << std::endl;

              std::cout << (Mdt - M0)/dMtStep << " & \\;(M_t - " << KPV.Mt() << ") & + ";
              std::cout << (MdH - M0)/dMHStep << " & \\;(M_H - " << KPV.MH() << ") & + ";
              std::cout << (MdW - M0)/dMWStep << " & \\;(M_W - " << KPV.MW() << ") & + ";
              std::cout << M0 << " \\non \\\\" << std::endl;      


              long double Y0  = it->second[0]->y(apow,aspow);
              long double Ydt = it->second[1]->y(apow,aspow);
              long double YdH = it->second[2]->y(apow,aspow);
              long double YdW = it->second[3]->y(apow,aspow);

              std::cout << "Interpolation for Y" << it->first.first 
                        << "["<< apow << "," << aspow <<"] at  scale \\mu = " 
                        << sqrt(it->first.second) << std::endl;

              std::cout << (Ydt - Y0)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
              std::cout << (YdH - Y0)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
              std::cout << (YdW - Y0)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
              std::cout << Y0 << std::endl;      

              // std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
              // std::cout << "& X_" << it->first.first 
              //           << "^{"<< apow << "," << aspow <<"}(\\mu = " 
              //           << sqrt(it->first.second) << ") = " << std::endl;

              // std::cout << (Ydt - Y0)/dMtStep << " (M_t - " << KPV.Mt() << ") + ";
              // std::cout << (YdH - Y0)/dMHStep << " (M_H - " << KPV.MH() << ") + ";
              // std::cout << (YdW - Y0)/dMWStep << " (M_W - " << KPV.MW() << ") + ";
              // std::cout << Y0 << " \\non \\\\" << std::endl;      

              
            }
      // // return(0);

      // // 
      // //            TT
      // // 
      // // Top quark mu=Mz
      // tt xMtt0      = tt(KPV  , KPV.MMZ());
      // tt xMttdivMt  = tt(divMt, KPV.MMZ());
      // tt xMttdivMH  = tt(divMH, KPV.MMZ());
      // tt xMttdivMW  = tt(divMW, KPV.MMZ());

      // long double tt1_10 = xMtt0.m10().real();
      // long double tt2dMt_10 = xMttdivMt.m10().real();
      // long double tt2dMH_10 = xMttdivMH.m10().real();
      // long double tt2dMW_10 = xMttdivMW.m10().real();


      // std::cout << std::endl;
      // std::cout << " Top quark \\mu=MZ d(Mt) = " << dMtStep << ", d(MH) = " << dMHStep << ", d(MW) = " << dMWStep << ":" << std::endl;
      // std::cout << "           X(tt)_1,0 = ";
      // std::cout << (tt2dMt_10 - tt1_10)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (tt2dMH_10 - tt1_10)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (tt2dMW_10 - tt1_10)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << tt1_10 << std::endl;      

      // long double tt1_11 = xMtt0.m11().real();
      // long double tt2dMt_11 = xMttdivMt.m11().real();
      // long double tt2dMH_11 = xMttdivMH.m11().real();
      // long double tt2dMW_11 = xMttdivMW.m11().real();

      // std::cout << "           X(tt)_1,1 = ";
      // std::cout << (tt2dMt_11 - tt1_11)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (tt2dMH_11 - tt1_11)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (tt2dMW_11 - tt1_11)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << tt1_11 << std::endl;      

      // long double tt1_20 = xMtt0.m20().real();
      // long double tt2dMt_20 = xMttdivMt.m20().real();
      // long double tt2dMH_20 = xMttdivMH.m20().real();
      // long double tt2dMW_20 = xMttdivMW.m20().real();

      // std::cout << "           X(tt)_2,0 = ";
      // std::cout << (tt2dMt_20 - tt1_20)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (tt2dMH_20 - tt1_20)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (tt2dMW_20 - tt1_20)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << tt1_20 << std::endl;      

      // // Top quark mu=Mt
      // xMtt0      = tt(KPV  , KPV.MMt());
      // xMttdivMt  = tt(divMt, KPV.MMt());
      // xMttdivMH  = tt(divMH, KPV.MMt());
      // xMttdivMW  = tt(divMW, KPV.MMt());

      // tt1_10 = xMtt0.m10().real();
      // tt2dMt_10 = xMttdivMt.m10().real();
      // tt2dMH_10 = xMttdivMH.m10().real();
      // tt2dMW_10 = xMttdivMW.m10().real();

      // std::cout << std::endl;
      // std::cout << " Top quark \\mu=Mt d(Mt) = " << dMtStep << ", d(MH) = " << dMHStep << ", d(MW) = " << dMWStep << ":" << std::endl;
      // std::cout << "           X(tt)_1,0 = ";
      // std::cout << (tt2dMt_10 - tt1_10)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (tt2dMH_10 - tt1_10)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (tt2dMW_10 - tt1_10)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << tt1_10 << std::endl;      

      // tt1_11 = xMtt0.m11().real();
      // tt2dMt_11 = xMttdivMt.m11().real();
      // tt2dMH_11 = xMttdivMH.m11().real();
      // tt2dMW_11 = xMttdivMW.m11().real();

      // std::cout << "           X(tt)_1,1 = ";
      // std::cout << (tt2dMt_11 - tt1_11)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (tt2dMH_11 - tt1_11)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (tt2dMW_11 - tt1_11)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << tt1_11 << std::endl;      

      // tt1_20 = xMtt0.m20().real();
      // tt2dMt_20 = xMttdivMt.m20().real();
      // tt2dMH_20 = xMttdivMH.m20().real();
      // tt2dMW_20 = xMttdivMW.m20().real();

      // std::cout << "           X(tt)_2,0 = ";
      // std::cout << (tt2dMt_20 - tt1_20)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (tt2dMH_20 - tt1_20)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (tt2dMW_20 - tt1_20)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << tt1_20 << std::endl;      










      // // 
      // //            WW
      // // 
      // // W boson mu=Mz
      // WW xMWW0      = WW(KPV  , KPV.MMZ());
      // WW xMWWdivMt  = WW(divMt, KPV.MMZ());
      // WW xMWWdivMH  = WW(divMH, KPV.MMZ());
      // WW xMWWdivMW  = WW(divMW, KPV.MMZ());

      // long double WW1_10 = xMWW0.m10().real();
      // long double WW2dMt_10 = xMWWdivMt.m10().real();
      // long double WW2dMH_10 = xMWWdivMH.m10().real();
      // long double WW2dMW_10 = xMWWdivMW.m10().real();


      // std::cout << std::endl;
      // std::cout << " W boson \\mu=MZ d(Mt) = " << dMtStep << ", d(MH) = " << dMHStep << ", d(MW) = " << dMWStep << ":" << std::endl;
      // std::cout << "           X(WW)_1,0 = ";
      // std::cout << (WW2dMt_10 - WW1_10)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (WW2dMH_10 - WW1_10)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (WW2dMW_10 - WW1_10)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << WW1_10 << std::endl;      

      // long double WW1_11 = xMWW0.m11().real();
      // long double WW2dMt_11 = xMWWdivMt.m11().real();
      // long double WW2dMH_11 = xMWWdivMH.m11().real();
      // long double WW2dMW_11 = xMWWdivMW.m11().real();

      // std::cout << "           X(WW)_1,1 = ";
      // std::cout << (WW2dMt_11 - WW1_11)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (WW2dMH_11 - WW1_11)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (WW2dMW_11 - WW1_11)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << WW1_11 << std::endl;      

      // long double WW1_20 = xMWW0.m20().real();
      // long double WW2dMt_20 = xMWWdivMt.m20().real();
      // long double WW2dMH_20 = xMWWdivMH.m20().real();
      // long double WW2dMW_20 = xMWWdivMW.m20().real();

      // std::cout << "           X(WW)_2,0 = ";
      // std::cout << (WW2dMt_20 - WW1_20)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (WW2dMH_20 - WW1_20)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (WW2dMW_20 - WW1_20)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << WW1_20 << std::endl;      

      // // W boson mu=Mt
      // xMWW0      = WW(KPV  , KPV.MMt());
      // xMWWdivMt  = WW(divMt, KPV.MMt());
      // xMWWdivMH  = WW(divMH, KPV.MMt());
      // xMWWdivMW  = WW(divMW, KPV.MMt());

      // WW1_10 = xMWW0.m10().real();
      // WW2dMt_10 = xMWWdivMt.m10().real();
      // WW2dMH_10 = xMWWdivMH.m10().real();
      // WW2dMW_10 = xMWWdivMW.m10().real();

      // std::cout << std::endl;
      // std::cout << " W boson \\mu=Mt d(Mt) = " << dMtStep << ", d(MH) = " << dMHStep << ", d(MW) = " << dMWStep << ":" << std::endl;
      // std::cout << "           X(WW)_1,0 = ";
      // std::cout << (WW2dMt_10 - WW1_10)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (WW2dMH_10 - WW1_10)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (WW2dMW_10 - WW1_10)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << WW1_10 << std::endl;      

      // WW1_11 = xMWW0.m11().real();
      // WW2dMt_11 = xMWWdivMt.m11().real();
      // WW2dMH_11 = xMWWdivMH.m11().real();
      // WW2dMW_11 = xMWWdivMW.m11().real();

      // std::cout << "           X(WW)_1,1 = ";
      // std::cout << (WW2dMt_11 - WW1_11)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (WW2dMH_11 - WW1_11)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (WW2dMW_11 - WW1_11)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << WW1_11 << std::endl;      

      // WW1_20 = xMWW0.m20().real();
      // WW2dMt_20 = xMWWdivMt.m20().real();
      // WW2dMH_20 = xMWWdivMH.m20().real();
      // WW2dMW_20 = xMWWdivMW.m20().real();

      // std::cout << "           X(WW)_2,0 = ";
      // std::cout << (WW2dMt_20 - WW1_20)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (WW2dMH_20 - WW1_20)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (WW2dMW_20 - WW1_20)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << WW1_20 << std::endl;      








      
      // // 
      // //            ZZ
      // // 
      // // Z boson mu=Mz
      // ZZ xMZZ0      = ZZ(KPV  , KPV.MMZ());
      // ZZ xMZZdivMt  = ZZ(divMt, KPV.MMZ());
      // ZZ xMZZdivMH  = ZZ(divMH, KPV.MMZ());
      // ZZ xMZZdivMW  = ZZ(divMW, KPV.MMZ());

      // long double ZZ1_10 = xMZZ0.m10().real();
      // long double ZZ2dMt_10 = xMZZdivMt.m10().real();
      // long double ZZ2dMH_10 = xMZZdivMH.m10().real();
      // long double ZZ2dMW_10 = xMZZdivMW.m10().real();


      // std::cout << std::endl;
      // std::cout << " Z boson \\mu=MZ d(Mt) = " << dMtStep << ", d(MH) = " << dMHStep << ", d(MW) = " << dMWStep << ":" << std::endl;
      // std::cout << "           X(ZZ)_1,0 = ";
      // std::cout << (ZZ2dMt_10 - ZZ1_10)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (ZZ2dMH_10 - ZZ1_10)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (ZZ2dMW_10 - ZZ1_10)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << ZZ1_10 << std::endl;      

      // long double ZZ1_11 = xMZZ0.m11().real();
      // long double ZZ2dMt_11 = xMZZdivMt.m11().real();
      // long double ZZ2dMH_11 = xMZZdivMH.m11().real();
      // long double ZZ2dMW_11 = xMZZdivMW.m11().real();

      // std::cout << "           X(ZZ)_1,1 = ";
      // std::cout << (ZZ2dMt_11 - ZZ1_11)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (ZZ2dMH_11 - ZZ1_11)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (ZZ2dMW_11 - ZZ1_11)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << ZZ1_11 << std::endl;      

      // long double ZZ1_20 = xMZZ0.m20().real();
      // long double ZZ2dMt_20 = xMZZdivMt.m20().real();
      // long double ZZ2dMH_20 = xMZZdivMH.m20().real();
      // long double ZZ2dMW_20 = xMZZdivMW.m20().real();

      // std::cout << "           X(ZZ)_2,0 = ";
      // std::cout << (ZZ2dMt_20 - ZZ1_20)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (ZZ2dMH_20 - ZZ1_20)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (ZZ2dMW_20 - ZZ1_20)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << ZZ1_20 << std::endl;      

      // // Z boson mu=Mt
      // xMZZ0      = ZZ(KPV  , KPV.MMt());
      // xMZZdivMt  = ZZ(divMt, KPV.MMt());
      // xMZZdivMH  = ZZ(divMH, KPV.MMt());
      // xMZZdivMW  = ZZ(divMW, KPV.MMt());

      // ZZ1_10 = xMZZ0.m10().real();
      // ZZ2dMt_10 = xMZZdivMt.m10().real();
      // ZZ2dMH_10 = xMZZdivMH.m10().real();
      // ZZ2dMW_10 = xMZZdivMW.m10().real();

      // std::cout << std::endl;
      // std::cout << " Z boson \\mu=Mt d(Mt) = " << dMtStep << ", d(MH) = " << dMHStep << ", d(MW) = " << dMWStep << ":" << std::endl;
      // std::cout << "           X(ZZ)_1,0 = ";
      // std::cout << (ZZ2dMt_10 - ZZ1_10)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (ZZ2dMH_10 - ZZ1_10)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (ZZ2dMW_10 - ZZ1_10)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << ZZ1_10 << std::endl;      

      // ZZ1_11 = xMZZ0.m11().real();
      // ZZ2dMt_11 = xMZZdivMt.m11().real();
      // ZZ2dMH_11 = xMZZdivMH.m11().real();
      // ZZ2dMW_11 = xMZZdivMW.m11().real();

      // std::cout << "           X(ZZ)_1,1 = ";
      // std::cout << (ZZ2dMt_11 - ZZ1_11)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (ZZ2dMH_11 - ZZ1_11)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (ZZ2dMW_11 - ZZ1_11)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << ZZ1_11 << std::endl;      

      // ZZ1_20 = xMZZ0.m20().real();
      // ZZ2dMt_20 = xMZZdivMt.m20().real();
      // ZZ2dMH_20 = xMZZdivMH.m20().real();
      // ZZ2dMW_20 = xMZZdivMW.m20().real();

      // std::cout << "           X(ZZ)_2,0 = ";
      // std::cout << (ZZ2dMt_20 - ZZ1_20)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (ZZ2dMH_20 - ZZ1_20)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (ZZ2dMW_20 - ZZ1_20)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << ZZ1_20 << std::endl;      



      // // 
      // //            HH
      // // 
      // // H boson mu=Mz
      // HH xMHH0      = HH(KPV  , KPV.MMZ());
      // HH xMHHdivMt  = HH(divMt, KPV.MMZ());
      // HH xMHHdivMH  = HH(divMH, KPV.MMZ());
      // HH xMHHdivMW  = HH(divMW, KPV.MMZ());

      // long double HH1_10 = xMHH0.m10().real();
      // long double HH2dMt_10 = xMHHdivMt.m10().real();
      // long double HH2dMH_10 = xMHHdivMH.m10().real();
      // long double HH2dMW_10 = xMHHdivMW.m10().real();


      // std::cout << std::endl;
      // std::cout << " Higgs \\mu=MZ d(Mt) = " << dMtStep << ", d(MH) = " << dMHStep << ", d(MW) = " << dMWStep << ":" << std::endl;
      // std::cout << "           X(HH)_1,0 = ";
      // std::cout << (HH2dMt_10 - HH1_10)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (HH2dMH_10 - HH1_10)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (HH2dMW_10 - HH1_10)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << HH1_10 << std::endl;      

      // long double HH1_11 = xMHH0.m11().real();
      // long double HH2dMt_11 = xMHHdivMt.m11().real();
      // long double HH2dMH_11 = xMHHdivMH.m11().real();
      // long double HH2dMW_11 = xMHHdivMW.m11().real();

      // std::cout << "           X(HH)_1,1 = ";
      // std::cout << (HH2dMt_11 - HH1_11)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (HH2dMH_11 - HH1_11)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (HH2dMW_11 - HH1_11)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << HH1_11 << std::endl;      

      // long double HH1_20 = xMHH0.m20().real();
      // long double HH2dMt_20 = xMHHdivMt.m20().real();
      // long double HH2dMH_20 = xMHHdivMH.m20().real();
      // long double HH2dMW_20 = xMHHdivMW.m20().real();

      // std::cout << "           X(HH)_2,0 = ";
      // std::cout << (HH2dMt_20 - HH1_20)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (HH2dMH_20 - HH1_20)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (HH2dMW_20 - HH1_20)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << HH1_20 << std::endl;      

      // // H boson mu=Mt
      // xMHH0      = HH(KPV  , KPV.MMt());
      // xMHHdivMt  = HH(divMt, KPV.MMt());
      // xMHHdivMH  = HH(divMH, KPV.MMt());
      // xMHHdivMW  = HH(divMW, KPV.MMt());

      // HH1_10 = xMHH0.m10().real();
      // HH2dMt_10 = xMHHdivMt.m10().real();
      // HH2dMH_10 = xMHHdivMH.m10().real();
      // HH2dMW_10 = xMHHdivMW.m10().real();

      // std::cout << std::endl;
      // std::cout << " Higgs \\mu=Mt d(Mt) = " << dMtStep << ", d(MH) = " << dMHStep << ", d(MW) = " << dMWStep << ":" << std::endl;
      // std::cout << "           X[HH,1,0]:= ";
      // std::cout << (HH2dMt_10 - HH1_10)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (HH2dMH_10 - HH1_10)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (HH2dMW_10 - HH1_10)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << HH1_10 << std::endl;      

      // HH1_11 = xMHH0.m11().real();
      // HH2dMt_11 = xMHHdivMt.m11().real();
      // HH2dMH_11 = xMHHdivMH.m11().real();
      // HH2dMW_11 = xMHHdivMW.m11().real();

      // std::cout << "           X[HH,1,1]:= ";
      // std::cout << (HH2dMt_11 - HH1_11)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (HH2dMH_11 - HH1_11)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (HH2dMW_11 - HH1_11)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << HH1_11 << std::endl;      

      // HH1_20 = xMHH0.m20().real();
      // HH2dMt_20 = xMHHdivMt.m20().real();
      // HH2dMH_20 = xMHHdivMH.m20().real();
      // HH2dMW_20 = xMHHdivMW.m20().real();

      // std::cout << "           X[HH,2,0]:= ";
      // std::cout << (HH2dMt_20 - HH1_20)/dMtStep << " * (Mt - " << KPV.Mt() << ") + ";
      // std::cout << (HH2dMH_20 - HH1_20)/dMHStep << " * (MH - " << KPV.MH() << ") + ";
      // std::cout << (HH2dMW_20 - HH1_20)/dMWStep << " * (MW - " << KPV.MW() << ") + ";
      // std::cout << HH1_20 << std::endl;      

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

