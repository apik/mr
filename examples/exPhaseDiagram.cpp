#include <iostream>
#include <cmath>
#include <fstream>
#include "mr.hpp"
#include "tools.hpp"
#include "boost/tuple/tuple.hpp"

using boost::tie;


bool pointsMH(std::string fname, long double scalePow, long double start, long double end, long double step = 1., long double inMtop = pdg2014::Mt, long double asMZ = pdg2014::asMZ)
{
  std::ofstream fout_bet(fname.c_str());
  
  if ( ((start < end) && step < 0) || ((start > end) && step > 0))
    return false;
  
  for (double MH = start; MH <= end; MH += step)
    {
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, MH, inMtop);
      
      double cMt; int state;
      // cMt = critMt(oi, pow(10,scalePow));
      tie(cMt, state)= critMt0<Mt_Beta0> (oi, pow(10,scalePow), asMZ);
      
      fout_bet << MH << " " << cMt << " " << state << std::endl;
      std::cout << "Crit Mt = " << cMt << " Mh = " << MH << " mu = 10^" << scalePow << std::endl;
      // std::cout << "New Mt" << critMt(oi, pow(10,scalePow)) << " for MH = " << MH << std::endl; 
    }
  return true;
}


bool pointsMt(std::string fname, long double scalePow, long double start, long double end, long double step = -1., long double inMH = pdg2014::MH, long double asMZ = pdg2014::asMZ)
{ 
  std::ofstream fout_bet(fname.c_str());
  
  if ( ((start < end) && step < 0) || ((start > end) && step > 0))
    return false;
  
  for (double Mt = start; Mt >= end; Mt += step)
    {
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, inMH, Mt);
      
      double cMH; int state;
      tie(cMH, state)= critMH0<MH_Beta0> (oi, pow(10,scalePow), asMZ);
      
      // cMH = critMH(oi, pow(10,scalePow));
      fout_bet << cMH << " " << Mt << " " << state << std::endl;
      
      std::cout << "Crit MH = " << cMH << " Mt = " << Mt << " mu = 10^" << scalePow << std::endl;
      // std::cout << "New MH " << cMH << " for Mt = " << Mt << std::endl; 
    }
  return 0;
}


bool pointsMHlam0(std::string fname, long double start, long double end, long double step = 1., long double asMZ = pdg2014::asMZ)
{
  std::ofstream fout_bet(fname.c_str());
  
  if ( ((start < end) && step < 0) || ((start > end) && step > 0))
    return false;
  
  for (double MH = start; MH <= end; MH += step)
    {
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, MH, pdg2014::Mt);
      
      double csc; double cPL; int state;
      // cMt = critMt(oi, pow(10,scalePow));
      
      fout_bet << MH << " ";
      for (size_t  scPow = 7; scPow <= 9; scPow++)
        {
          tie(csc, state) = critMt0<Mt_Lambda0> (oi, pow(10,scPow), asMZ);
          
          fout_bet << csc << " ";
        }

      tie(cPL, state) = critMt0<Mt_Lambda0> (oi, pdg2014::Mpl, asMZ);
      fout_bet << cPL << std::endl;
    }
  
  return true;
}

int main (int argc, char *argv[])
{
  try
    {


      // Default log level is logERROR
      loglevel = logINFO;

      // #define instability      
#ifdef instability
      std::pair<long double,long double> metaStable_Instable(160, pow(10,8));
      

      
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, metaStable_Instable.first);
      
      std::ofstream finstable("instable.dat");
      
      for (double MH = 120.; MH <= 132.; MH += 2.)
        {
          
          oi.setMH(MH);
          metaStable_Instable = critMt_scaleNotFixed<Mt_Instability>(oi, metaStable_Instable.second, pdg2014::asMZ);
          
          std::cout << "\n\n\n\nCritical Top mass    at Mh = " << oi.MH() << " is " 
                    << metaStable_Instable.first   << std::endl; 
          
          
          finstable << MH << " " 
                    << metaStable_Instable.first << " " 
                    << metaStable_Instable.second << std::endl;
        }      

      
#else
     
      // Start point
      std::pair<long double,long double> metaStable_Stable_m(pdg2014::Mt, pdg2014::Mpl);
      std::pair<long double,long double> metaStable_Stable_c(pdg2014::Mt, pdg2014::Mpl);
      std::pair<long double,long double> metaStable_Stable_p(pdg2014::Mt, pdg2014::Mpl);
      
      
      
      OSinput oi(pdg2014::Mb, pdg2014::MW, pdg2014::MZ, pdg2014::MH, metaStable_Stable_c.first);
      
      std::ofstream fstable("stable.dat");

      for (double MH = 120.; MH <= 140.; MH += 2.)
        {
      
          oi.setMH(MH);
          metaStable_Stable_m = critMt_scaleNotFixed<Mt_Stability>(oi, metaStable_Stable_c.second, pdg2014::asMZ - 0.0006);
          metaStable_Stable_c = critMt_scaleNotFixed<Mt_Stability>(oi, metaStable_Stable_c.second, pdg2014::asMZ);
          metaStable_Stable_p = critMt_scaleNotFixed<Mt_Stability>(oi, metaStable_Stable_c.second, pdg2014::asMZ + 0.0006);

          std::cout << "\n\n\n\nCritical Top mass    at Mh = " << oi.MH() << " is " 
                    << metaStable_Stable_c.first   << " [+] " 
                    << (metaStable_Stable_p.first - metaStable_Stable_c.first) << "[-] " 
                    << (metaStable_Stable_m.first - metaStable_Stable_c.first) << std::endl; 


          fstable << MH << " " 
                  << metaStable_Stable_c.first << " " << metaStable_Stable_c.second << " "
                  << metaStable_Stable_p.first << " " << metaStable_Stable_p.second << " " 
                  << metaStable_Stable_m.first << " " << metaStable_Stable_m.second << std::endl;
        }      
#endif
      

      // MU=10^17
      // pointsMt("plotZeroLamBet/b17tm.dat", 17., 180, 166, -0.5, pdg2014::MH, pdg2014::asMZ - 0.0006);
      // pointsMt("plotZeroLamBet/b17tc.dat", 17., 180, 166, -0.5, pdg2014::MH, pdg2014::asMZ);
      // pointsMt("plotZeroLamBet/b17tp.dat", 17., 180, 166, -0.5, pdg2014::MH, pdg2014::asMZ + 0.0006);

      // pointsMH("plotZeroLamBet/b17hm.dat", 17., 122, 140, 0.5, 166, pdg2014::asMZ - 0.0006);
      // pointsMH("plotZeroLamBet/b17hc.dat", 17., 122, 140, 0.5, 166, pdg2014::asMZ);
      // pointsMH("plotZeroLamBet/b17hp.dat", 17., 122, 140, 0.5, 166, pdg2014::asMZ + 0.0006);



      // MU=10^18
      // pointsMt("plotZeroLamBet/b18tm.dat", 18., 180, 169, -0.5, pdg2014::MH, pdg2014::asMZ - 0.0006);
      // pointsMt("plotZeroLamBet/b18tc.dat", 18., 180, 169, -0.5, pdg2014::MH, pdg2014::asMZ);
      // pointsMt("plotZeroLamBet/b18tp.dat", 18., 180, 169, -0.5, pdg2014::MH, pdg2014::asMZ + 0.0006);

      // pointsMH("plotZeroLamBet/b18hm.dat", 18., 127, 140, 0.5, 166, pdg2014::asMZ - 0.0006);
      // pointsMH("plotZeroLamBet/b18hc.dat", 18., 127, 140, 0.5, 166, pdg2014::asMZ);
      // pointsMH("plotZeroLamBet/b18hp.dat", 18., 127, 140, 0.5, 166, pdg2014::asMZ + 0.0006);



      // MU=10^19
      // pointsMt("plotZeroLamBet/b19tm.dat", 19., 180, 171, -0.5, pdg2014::MH, pdg2014::asMZ - 0.0006);
      // pointsMt("plotZeroLamBet/b19tc.dat", 19., 180, 171, -0.5, pdg2014::MH, pdg2014::asMZ);
      // pointsMt("plotZeroLamBet/b19tp.dat", 19., 180, 171, -0.5, pdg2014::MH, pdg2014::asMZ + 0.0006);

      // pointsMH("plotZeroLamBet/b19hm.dat", 19., 131, 140, 0.5, 166, pdg2014::asMZ - 0.0006);
      // pointsMH("plotZeroLamBet/b19hc.dat", 19., 131, 140, 0.5, 166, pdg2014::asMZ);
      // pointsMH("plotZeroLamBet/b19hp.dat", 19., 131, 140, 0.5, 166, pdg2014::asMZ + 0.0006);




      // MU=M_PL
      // pointsMt("plotZeroLamBet/bPLtm.dat", log10(pdg2014::Mpl), 180, 172, -0.5, pdg2014::MH, pdg2014::asMZ - 0.0006);
      // pointsMt("plotZeroLamBet/bPLtc.dat", log10(pdg2014::Mpl), 180, 172, -0.5, pdg2014::MH, pdg2014::asMZ);
      // pointsMt("plotZeroLamBet/bPLtp.dat", log10(pdg2014::Mpl), 180, 172, -0.5, pdg2014::MH, pdg2014::asMZ + 0.0006);

      // pointsMH("plotZeroLamBet/bPLhm.dat", log10(pdg2014::Mpl), 131, 140, 0.5, 166, pdg2014::asMZ - 0.0006);
      // pointsMH("plotZeroLamBet/bPLhc.dat", log10(pdg2014::Mpl), 131, 140, 0.5, 166, pdg2014::asMZ);
      // pointsMH("plotZeroLamBet/bPLhp.dat", log10(pdg2014::Mpl), 131, 140, 0.5, 166, pdg2014::asMZ + 0.0006);




      

      // Zeros of Lambda
      // pointsMHlam0("plotZeroLamBet/lam0m.dat", 120, 140, 2, pdg2014::asMZ - 0.0006);
      pointsMHlam0("plotZeroLamBet/lam789.dat", 120, 140, 2, pdg2014::asMZ);
      // pointsMHlam0("plotZeroLamBet/lam0p.dat", 120, 140, 2, pdg2014::asMZ + 0.0006);



    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

