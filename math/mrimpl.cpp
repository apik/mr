#include <iostream>
#include <map>
#include <mathlink.h>
#include "bb.hpp"
#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp"
#include "dr.hpp"


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

const dr<MS> & get_drbar(const MSinput& mi, long double mu2)
{
  static std::map< std::pair<MSinput, long double>, dr<MS> > directory;
  std::map< std::pair<MSinput, long double>, dr<MS> >::iterator i = directory.find(std::make_pair(mi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(mi,mu2), dr<MS>(mi,mu2))).first->second;
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

const bb & get_bb(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, bb > directory;
  std::map< std::pair<OSinput, long double>, bb >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), bb(oi,mu2))).first->second;
}

typedef std::vector<std::pair<size_t,size_t> > Vpow;


void MW(long double gp, long double g, long double gs, long double yt, long double lam, long double mu0, long double mu, int L) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromConsts(mu2, mu0, lam, 0, yt, g, gp);
  
  WW<MS> WWm = get_WWbar(mi, mu2);

  long double aEW  = mi.alpha()/4./Pi;
  long double aQCD = pow(gs/4./Pi,2);

  long double dMW2 = 0;
  if (L>0) dMW2 += aEW*WWm.x10();
  if (L>1) dMW2 += aEW*aQCD*WWm.x11() 
	  	  +aEW*aEW*WWm.x20();

  long double MWpole = mi.mW()*sqrt(1 + dMW2);
  MLPutReal128(stdlink, MWpole);
  
}

void MZ(long double gp, long double g, long double gs, long double yt, long double lam, long double mu0, long double mu, int L) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromConsts(mu2, mu0, lam, 0, yt, g, gp);
  
  ZZ<MS> ZZm = get_ZZbar(mi, mu2);

  long double aEW  = mi.alpha()/4./Pi;
  long double aQCD = pow(gs/4./Pi,2);

  long double dMZ2 = 0;
  if (L>0) dMZ2 += aEW*ZZm.x10();
  if (L>1) dMZ2 += aEW*aQCD*ZZm.x11() 
	  	  +aEW*aEW*ZZm.x20();

  long double MZpole = mi.mZ()*sqrt(1 + dMZ2);
  MLPutReal128(stdlink, MZpole);
  
}

void MH(long double gp, long double g, long double gs, long double yt, long double lam, long double mu0, long double mu, int L) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromConsts(mu2, mu0, lam, 0, yt, g, gp);
  
  HH<MS> HHm = get_HHbar(mi, mu2);

  long double aEW  = mi.alpha()/4./Pi;
  long double aQCD = pow(gs/4./Pi,2);

  long double dMH2 = 0;
  if (L>0) dMH2 += aEW*HHm.x10();
  if (L>1) dMH2 += aEW*aQCD*HHm.x11() 
	  	  +aEW*aEW*HHm.x20();

  long double MHpole = mi.mH()*sqrt(1 + dMH2);
  MLPutReal128(stdlink, MHpole);
  
}
void MT(long double gp, long double g, long double gs, long double yt, long double lam, long double mu0, long double mu, int L) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromConsts(mu2, mu0, lam, 0, yt, g, gp);
  
  tt<MS> ttm = get_ttbar(mi, mu2);

  long double aEW  = mi.alpha()/4./Pi;
  long double aQCD = pow(gs/4./Pi,2);

  long double dMT2 = 0;
  if (L>0) dMT2 += aEW*ttm.x10() + aQCD*ttm.x01();
  if (L>1) dMT2 += aEW*aQCD*ttm.x11() 
	  	  +aEW*aEW*ttm.x20()
	          +aQCD*aQCD*ttm.x02();
  if (L>2) dMT2 +=  aQCD*aQCD*aQCD*ttm.x03();

  long double MTpole = mi.mt()*(1 + dMT2);

  MLPutReal128(stdlink, MTpole);
  
}
void GF(long double gp, long double g, long double gs, long double yt, long double lam, long double mu0, long double mu, int L) 
{
  long double mu2 = pow(mu,2);
  MSinput mi = MSinput::fromConsts(mu2, mu0, lam, 0, yt, g, gp);
  
  dr<MS> drm = get_drbar(mi, mu2);

  long double aEW  = mi.alpha()/4./Pi;
  long double aQCD = pow(gs/4./Pi,2);

  long double dGF = 0;
  if (L>0) dGF += aEW*drm.dr10();
  if (L>1) dGF += aEW*aQCD*drm.dr11() 
	  	  +aEW*aEW*drm.dr20();

  long double GF = (1 + dGF)/sqrt(2.0)/pow(mi.vev(),2);
  MLPutReal128(stdlink, GF);
  
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
        MLPutFunction(stdlink, "yt", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, ttm.y(apow, aspow, nL, nH));
      }
}

void Xb(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  bb bbm = get_bb(oi, pow(mu,2));

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
        MLPutFunction(stdlink, "yb", 2);
        MLPutInteger(stdlink, apow);
        MLPutInteger(stdlink, aspow);
        MLPutReal128(stdlink, bbm.y(apow, aspow, nL, nH));
      }
}


// Pure QCD corrections

void XbQCD(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nl, int nh) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  bb bbm = get_bb(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 3);

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

}

void XtQCD(long double mb, long double mW, long double mZ, long double mH, long double mt, long double mu, int nl, int nh) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  tt<OS> ttm = get_tt(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 3);

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

}

int main(int argc, char* argv[]) 
{
  return MLMain(argc, argv);
}
