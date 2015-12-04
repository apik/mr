#include <iostream>
#include <map>
#include <mathlink.h>
#include "bb.hpp"
#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp"
#include "dr.hpp"
#include "alphaGF.hpp"
#include "smRGE.hpp"
#include "alphas.hpp"

using namespace mr;

typedef CouplingsSM<3,3,3,3,3,3,3> SM_3;
typedef CouplingsSM<2,2,2,2,2,2,2> SM_2;
typedef CouplingsSM<1,1,1,1,1,1,1> SM_1;

const WW<MS> & get_WWbar(const MSinput& mi, long double mu2)
{
  static std::map< std::pair<MSinput, long double>, WW<MS> > directory;
  std::map< std::pair<MSinput, long double>, WW<MS> >::iterator i = directory.find(std::make_pair(mi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(mi,mu2), WW<MS>(mi,mu2))).first->second;
}

const ZZ<MS> & get_ZZbar(const MSinput& mi, long double mu2)
{
  static std::map< std::pair<MSinput, long double>, ZZ<MS> > directory;
  std::map< std::pair<MSinput, long double>, ZZ<MS> >::iterator i = directory.find(std::make_pair(mi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(mi,mu2), ZZ<MS>(mi,mu2))).first->second;
}

const HH<MS> & get_HHbar(const MSinput& mi, long double mu2)
{
  static std::map< std::pair<MSinput, long double>, HH<MS> > directory;
  std::map< std::pair<MSinput, long double>, HH<MS> >::iterator i = directory.find(std::make_pair(mi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(mi,mu2), HH<MS>(mi,mu2))).first->second;
}

const tt<MS> & get_ttbar(const MSinput& mi, long double mu2)
{
  static std::map< std::pair<MSinput, long double>, tt<MS> > directory;
  std::map< std::pair<MSinput, long double>, tt<MS> >::iterator i = directory.find(std::make_pair(mi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(mi,mu2), tt<MS>(mi,mu2))).first->second;
}

const bb<MS> & get_bbbar(const MSinput& mi, long double mu2)
{
  static std::map< std::pair<MSinput, long double>, bb<MS> > directory;
  std::map< std::pair<MSinput, long double>, bb<MS> >::iterator i = directory.find(std::make_pair(mi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(mi,mu2), bb<MS>(mi,mu2))).first->second;
}
const dr<MS> & get_drbar(const MSinput& mi, long double mu2)
{
  static std::map< std::pair<MSinput, long double>, dr<MS> > directory;
  std::map< std::pair<MSinput, long double>, dr<MS> >::iterator i = directory.find(std::make_pair(mi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(mi,mu2), dr<MS>(mi,mu2))).first->second;
}

const dr<OS> & get_dr(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, dr<OS> > directory;
  std::map< std::pair<OSinput, long double>, dr<OS> >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), dr<OS>(oi,mu2))).first->second;
}

const WW<OS> & get_WW(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, WW<OS> > directory;
  std::map< std::pair<OSinput, long double>, WW<OS> >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), WW<OS>(oi,mu2))).first->second;
}

const ZZ<OS> & get_ZZ(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, ZZ<OS> > directory;
  std::map< std::pair<OSinput, long double>, ZZ<OS> >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), ZZ<OS>(oi,mu2))).first->second;
}

const HH<OS> & get_HH(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, HH<OS> > directory;
  std::map< std::pair<OSinput, long double>, HH<OS> >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), HH<OS>(oi,mu2))).first->second;
}

const tt<OS> & get_tt(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, tt<OS> > directory;
  std::map< std::pair<OSinput, long double>, tt<OS> >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), tt<OS>(oi,mu2))).first->second;
}

const bb<OS> & get_bb(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, bb<OS> > directory;
  std::map< std::pair<OSinput, long double>, bb<OS> >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), bb<OS>(oi,mu2))).first->second;
}

const alphaGF & get_aGF(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, alphaGF > directory;
  std::map< std::pair<OSinput, long double>, alphaGF >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), alphaGF(oi,mu2))).first->second;
}

typedef std::vector<std::pair<size_t,size_t> > Vpow;

void XMMW(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double mu) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromCouplings(g1, g2, yb, yt, lam, mu0, mu);
  
  WW<MS> WWm = get_WWbar(mi, mu2);


  MLPutFunction(stdlink, "List", 3);


  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
	{
		MLPutFunction(stdlink, "Rule", 2);
		MLPutFunction(stdlink, "xMMW", 2);
		MLPutInteger(stdlink, apow);
		MLPutInteger(stdlink, aspow);
  		MLPutReal128(stdlink, WWm.x(apow,aspow));
	}
 }

void XMMZ(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double mu) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromCouplings(g1, g2, yb, yt, lam, mu0, mu);
  
  ZZ<MS> ZZm = get_ZZbar(mi, mu2);


  MLPutFunction(stdlink, "List", 3);


  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
	{
		MLPutFunction(stdlink, "Rule", 2);
		MLPutFunction(stdlink, "xMMZ", 2);
		MLPutInteger(stdlink, apow);
		MLPutInteger(stdlink, aspow);
  		MLPutReal128(stdlink, ZZm.x(apow,aspow));
	}
}
void XMMH(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double mu) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromCouplings(g1, g2, yb, yt, lam, mu0, mu);
  
  HH<MS> HHm = get_HHbar(mi, mu2);


  MLPutFunction(stdlink, "List", 3);


  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
	{
		MLPutFunction(stdlink, "Rule", 2);
		MLPutFunction(stdlink, "xMMH", 2);
		MLPutInteger(stdlink, apow);
		MLPutInteger(stdlink, aspow);
  		MLPutReal128(stdlink, HHm.x(apow,aspow));
	}
}


void XdRbar(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double mu) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromCouplings(g1, g2, yb, yt, lam, mu0, mu);
  
  dr<MS> drm = get_drbar(mi, mu2);

  MLPutFunction(stdlink, "List", 3);

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xdRbar", 2);
  MLPutInteger(stdlink, 1);
  MLPutInteger(stdlink, 0);
  MLPutReal128(stdlink, drm.dr10());

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xdRbar", 2);
  MLPutInteger(stdlink, 1);
  MLPutInteger(stdlink, 1);
  MLPutReal128(stdlink, drm.dr11());

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xdRbar", 2);
  MLPutInteger(stdlink, 2);
  MLPutInteger(stdlink, 0);
  MLPutReal128(stdlink, drm.dr20());

}


void XMT(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double mu) 
{
  long double mu2 = pow(mu,2);

  MSinput mi = MSinput::fromCouplings(g1, g2, yb, yt, lam, mu0, mu);
  
  tt<MS> ttm = get_ttbar(mi, mu2);


  MLPutFunction(stdlink, "List", 3);


  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
	{
		MLPutFunction(stdlink, "Rule", 2);
		MLPutFunction(stdlink, "xMT", 2);
		MLPutInteger(stdlink, apow);
		MLPutInteger(stdlink, aspow);
  		MLPutReal128(stdlink, ttm.x(apow,aspow));
	}
 }


void XMTQCD(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double mu) 
{
  long double mu2 = pow(mu,2);

  MSinput mi = MSinput::fromCouplings(g1, g2, yb, yt, lam, mu0, mu);
  
  tt<MS> ttm = get_ttbar(mi, mu2);


  MLPutFunction(stdlink, "List", 3);


  for(size_t aspow = 1; aspow <=3; aspow++)
	{
		MLPutFunction(stdlink, "Rule", 2);
		MLPutFunction(stdlink, "xMT", 2);
		MLPutInteger(stdlink, 0);
		MLPutInteger(stdlink, aspow);
  		MLPutReal128(stdlink, ttm.x(0,aspow));
	}
 }

void XMB(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double mu) 
{
  long double mu2 = pow(mu,2);

  MSinput mi = MSinput::fromCouplings(g1, g2, yb, yt, lam, mu0, mu);
  
  bb<MS> bbm = get_bbbar(mi, mu2);


  MLPutFunction(stdlink, "List", 3);


  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
	{
		MLPutFunction(stdlink, "Rule", 2);
		MLPutFunction(stdlink, "xMB", 2);
		MLPutInteger(stdlink, apow);
		MLPutInteger(stdlink, aspow);
  		MLPutReal128(stdlink, bbm.x(apow,aspow));
	}
 }


void XMBQCD(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double mu) 
{
  long double mu2 = pow(mu,2);

  MSinput mi = MSinput::fromCouplings(g1, g2, yb, yt, lam, mu0, mu);
  
  bb<MS> bbm = get_bbbar(mi, mu2);


  MLPutFunction(stdlink, "List", 2);


  for(size_t aspow = 1; aspow <=2; aspow++)
	{
		MLPutFunction(stdlink, "Rule", 2);
		MLPutFunction(stdlink, "xMB", 2);
		MLPutInteger(stdlink, 0);
		MLPutInteger(stdlink, aspow);
  		MLPutReal128(stdlink, bbm.x(0,aspow));
	}
 }


void RunQCDnf6(long double oscale, long double asMZ, long double MZscale, int nL, long double mtpole, long double mtth)
{
	// nL - number of loops
	// mtth - top threshold
	// run to threshold - MT - 4 loop RGE
	long double asMT = run(asMZ, MZscale, mtth, 5, nL);
	// jump to nf = 6
	asMT = as5nf2as6nf(mtpole, mtth, asMT, 5, nL - 1); 
	MLPutReal128(stdlink,run(asMT, mtth, oscale, 6, nL));
		
}	
void RunQCD(long double oscale, long double as0, long double iscale, int nL, int nF)
{
	// nL - number of loops
	// mtth - top threshold
	// run to threshold - MT - 4 loop RGE
	long double aso = run(as0, iscale, oscale, nF, nL);
	// jump to nf = 6
	MLPutReal128(stdlink,aso);
		
}	

void RunSM(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double iscale, long double oscale) 
{
      Couplings<3,3,3,3,3,3,3,3,3> runSM(
                pow(g1/4./Pi,2),
                pow(g2/4./Pi,2),
                pow(gs/4./Pi,2),
                pow(yt/4./Pi,2),
                pow(yb/4./Pi,2), // yb?
                0, // ytau?
                lam/pow(4.*Pi,2),             // Lambda
		mu0, /* higgs mass parameter*/
		0, /* higgs vev */
                pow(iscale,2),
                3 // NG
                );

      SMCouplings runCoupling = runSM(pow(oscale,2));
      MLPutFunction(stdlink, "List", 8);
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[0])); //g1
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[1])); //g2 
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[2])); //gs
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[4])); //yb
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[3])); //yt
      MLPutReal128(stdlink, pow(4.*Pi,2)*runCoupling[6]); //lam
      MLPutReal128(stdlink, runCoupling[7]); //m ??
      MLPutReal128(stdlink, oscale); // out scale
		
}

void RunSMwithBetas(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double mu0, long double iscale, long double oscale) 
{
      Couplings<3,3,3,3,3,3,3,3,3> runSM(
                pow(g1/4./Pi,2),
                pow(g2/4./Pi,2),
                pow(gs/4./Pi,2),
                pow(yt/4./Pi,2),
                pow(yb/4./Pi,2), // yb?
                0, // ytau?
                lam/pow(4.*Pi,2),             // Lambda
		mu0, /* higgs mass parameter*/
		0, /* higgs vev */
                pow(iscale,2),
                3 // NG
                );

      std::pair<SMCouplings, SMCouplings> runCouplingAndBeta = runSM.AandB(pow(oscale,2));
      MLPutFunction(stdlink, "List", 8);
      MLPutFunction(stdlink, "List", 2);
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCouplingAndBeta.first[0]));  //g1
      MLPutReal128(stdlink, runCouplingAndBeta.second[0]);             //g1
      MLPutFunction(stdlink, "List", 2);
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCouplingAndBeta.first[1]));  //g2 
      MLPutReal128(stdlink, runCouplingAndBeta.second[1]);             //g2
      MLPutFunction(stdlink, "List", 2);
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCouplingAndBeta.first[2]));  //gs
      MLPutReal128(stdlink, runCouplingAndBeta.second[2]);             //gs
      MLPutFunction(stdlink, "List", 2);
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCouplingAndBeta.first[4]));  //yb
      MLPutReal128(stdlink, runCouplingAndBeta.second[4]);             //yb
      MLPutFunction(stdlink, "List", 2);
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCouplingAndBeta.first[3]));  //yt
      MLPutReal128(stdlink, runCouplingAndBeta.second[3]);             //yt
      MLPutFunction(stdlink, "List", 2);
      MLPutReal128(stdlink, pow(4.*Pi,2)*runCouplingAndBeta.first[6]); //lam
      MLPutReal128(stdlink, runCouplingAndBeta.second[6]);             //lam
      MLPutFunction(stdlink, "List", 2);
      MLPutReal128(stdlink, runCouplingAndBeta.first[7]);              //m0 ??
      MLPutReal128(stdlink, runCouplingAndBeta.second[7]);             //m0
      MLPutReal128(stdlink, oscale); // out scale
		
}

void RunSMcouplings(long double g1, long double g2, long double gs, long double yb, long double yt, long double lam, long double iscale, long double oscale, int loop) 
{

      SMCouplings runCoupling;

      switch(loop)
      {
	      case 1:
		{
		      SM_1 runSM1(
                		pow(g1/4./Pi,2),
                		pow(g2/4./Pi,2),
                		pow(gs/4./Pi,2),
                		pow(yt/4./Pi,2),
                		pow(yb/4./Pi,2), // yb?
                		0, // ytau?
                		lam/pow(4.*Pi,2),             // Lambda
                		pow(iscale,2),
                		3 // NG
                		);
		      runCoupling = runSM1(pow(oscale,2));
		}
		      break;
	      case 2:
		{
		      SM_2 runSM2(
                		pow(g1/4./Pi,2),
                		pow(g2/4./Pi,2),
                		pow(gs/4./Pi,2),
                		pow(yt/4./Pi,2),
                		pow(yb/4./Pi,2), // yb?
                		0, // ytau?
                		lam/pow(4.*Pi,2),             // Lambda
                		pow(iscale,2),
                		3 // NG
                		);
		      runCoupling = runSM2(pow(oscale,2));
		 }
		      break;
	      case 3:
		 {
		      SM_3 runSM3(
                		pow(g1/4./Pi,2),
                		pow(g2/4./Pi,2),
                		pow(gs/4./Pi,2),
                		pow(yt/4./Pi,2),
                		pow(yb/4./Pi,2), // yb?
                		0, // ytau?
                		lam/pow(4.*Pi,2),             // Lambda
                		pow(iscale,2),
                		3 // NG
                		);
		      runCoupling = runSM3(pow(oscale,2));
		 }
		      break;
	      default:
		 {
		      SM_3 runSM(
                		pow(g1/4./Pi,2),
                		pow(g2/4./Pi,2),
                		pow(gs/4./Pi,2),
                		pow(yt/4./Pi,2),
                		pow(yb/4./Pi,2), // yb?
                		0, // ytau?
                		lam/pow(4.*Pi,2),             // Lambda
                		pow(iscale,2),
                		3 // NG
                		);
		      runCoupling = runSM(pow(oscale,2));
		 }
	
	}	

      MLPutFunction(stdlink, "List", 7);
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[0])); //g1
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[1])); //g2 
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[2])); //gs
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[4])); //yb
      MLPutReal128(stdlink, 4.*Pi*sqrt(runCoupling[3])); //yt
      MLPutReal128(stdlink, pow(4.*Pi,2)*runCoupling[6]); //lam
      MLPutReal128(stdlink, oscale); // out scale
		
}



void XW(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  WW<OS> WWm = get_WW(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 6);

  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
      {
        // Mass
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "xW", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, WWm.x(apow, aspow, nL, nH));
        // Yukawa
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "yW", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, WWm.y(apow, aspow, nL, nH));
      }
}

void XZ(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  ZZ<OS> ZZm = get_ZZ(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 6);

  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
      {
        // Mass
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "xZ", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, ZZm.x(apow, aspow, nL, nH));
        // Yukawa
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "yZ", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, ZZm.y(apow, aspow, nL, nH));
      }
}

void XH(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  HH<OS> HHm = get_HH(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 6);

  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
      {
        // Mass
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "xH", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, HHm.x(apow, aspow, nL, nH));
        // Yukawa
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "yH", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, HHm.y(apow, aspow, nL, nH));
      }
}

void Xt(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  tt<OS> ttm = get_tt(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 6);

  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
      {
        // Mass
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "xt", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, ttm.x(apow, aspow, nL, nH));
        // Yukawa
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "yT", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, ttm.y(apow, aspow, nL, nH));
      }
}

void dROS(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  dr<OS> dros = get_dr(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 3);

        // dR
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "dr", 2);
        MLPutInteger(stdlink, 1);
        MLPutInteger(stdlink, 0);
        MLPutReal128(stdlink, dros.dr10());

        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "dr", 2);
        MLPutInteger(stdlink, 1);
        MLPutInteger(stdlink, 1);
        MLPutReal128(stdlink, dros.dr11());

        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "dr", 2);
        MLPutInteger(stdlink, 2);
        MLPutInteger(stdlink, 0);
        MLPutReal128(stdlink, dros.dr20());
}

void dalphaGF(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  alphaGF agf = get_aGF(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 3);

        // dalphaGF
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "daGF", 2);
        MLPutInteger(stdlink, 1);
        MLPutInteger(stdlink, 0);
        MLPutReal128(stdlink, agf.a10(nL,nH,1));

        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "daGF", 2);
        MLPutInteger(stdlink, 1);
        MLPutInteger(stdlink, 1);
        MLPutReal128(stdlink, agf.a11(nL,nH,1));

        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "daGF", 2);
        MLPutInteger(stdlink, 2);
        MLPutInteger(stdlink, 0);
        MLPutReal128(stdlink, agf.a20(nL,nH,1));

}

void Xb(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  bb<OS> bbm = get_bb(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 6);

  for(size_t apow = 1; apow <=2; apow++)
    for(size_t aspow = 0; aspow + apow <=2; aspow++)
      {
        // Mass
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "xb", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, bbm.x(apow, aspow, nL, nH));
        // Yukawa
        MLPutFunction(stdlink, "Rule", 2);
        MLPutFunction(stdlink, "yB", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, bbm.y(apow, aspow, nL, nH));
      }
}


// NB: probably wrong: Pure QCD corrections

void XbQCD(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nl, int nh) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  bb<OS> bbm = get_bb(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 3*2);

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xb", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 1);
  MLPutReal128(stdlink, bbm.x01(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xb", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 2);
  MLPutReal128(stdlink, bbm.x02(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xb", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 3);
  MLPutReal128(stdlink, bbm.x03(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "yB", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 1);
  MLPutReal128(stdlink, bbm.x01(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "yB", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 2);
  MLPutReal128(stdlink, bbm.x02(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "yB", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 3);
  MLPutReal128(stdlink, bbm.x03(nl,nh));
}

void XtQCD(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nl, int nh) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  tt<OS> ttm = get_tt(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 4*2);

// contrib to coupling to mmt/MMT
  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xt", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 1);
  MLPutReal128(stdlink, ttm.x01(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xt", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 2);
  MLPutReal128(stdlink, ttm.x02(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xt", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 3);
  MLPutReal128(stdlink, ttm.x03(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xt", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 4);
  MLPutReal128(stdlink, ttm.x04(nl,nh));
// contrib to coupling
  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "yT", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 1);
  MLPutReal128(stdlink, ttm.x01(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "yT", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 2);
  MLPutReal128(stdlink, ttm.x02(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "yT", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 3);
  MLPutReal128(stdlink, ttm.x03(nl,nh));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "yT", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 4);
  MLPutReal128(stdlink, ttm.x04(nl,nh));

}

int main(int argc, char* argv[]) 
{
  return MLMain(argc, argv);
}
