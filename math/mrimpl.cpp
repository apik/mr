#include <iostream>
#include <map>
#include <mathlink.h>
#include "bb.hpp"
#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp"



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

const tt & get_tt(const OSinput& oi, long double mu2)
{
  static std::map< std::pair<OSinput, long double>, tt > directory;
  std::map< std::pair<OSinput, long double>, tt >::iterator i = directory.find(std::make_pair(oi,mu2));
  if (i != directory.end())
    return i->second;
  else
    return directory.insert(std::make_pair(std::make_pair(oi,mu2), tt(oi,mu2))).first->second;
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
void XW(double mb, double mW, double mZ, double mH, double mt, double mu, int nL,int nH) 
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

void XZ(double mb, double mW, double mZ, double mH, double mt, double mu, int nL,int nH) 
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

void XH(double mb, double mW, double mZ, double mH, double mt, double mu, int nL,int nH) 
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

void Xt(double mb, double mW, double mZ, double mH, double mt, double mu, int nL,int nH) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  tt ttm = get_tt(oi, pow(mu,2));

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

void Xb(double mb, double mW, double mZ, double mH, double mt, double mu, int nL,int nH) 
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

// void XW(double mb, double mW, double mZ, double mH, double mt, double mu, int apow, int aspow, int nL,int nH) 
// {
//   OSinput oi(mb, mW, mZ, mH, mt);

//   WW<OS> wwm = get_WW(oi, pow(mu,2));

//   MLPutReal128(stdlink, wwm.x(apow, aspow, nL, nH));
// }

// void XZ(double mb, double mW, double mZ, double mH, double mt, double mu, int apow, int aspow, int nL,int nH) 
// {
//   OSinput oi(mb, mW, mZ, mH, mt);

//   ZZ<OS> zzm = get_ZZ(oi, pow(mu,2));

//   MLPutReal128(stdlink, zzm.x(apow, aspow, nL, nH));
// }

// void XH(double mb, double mW, double mZ, double mH, double mt, double mu, int apow, int aspow, int nL,int nH) 
// {
//   OSinput oi(mb, mW, mZ, mH, mt);

//   HH<OS> hhm = get_HH(oi, pow(mu,2));

//   MLPutReal128(stdlink, hhm.x(apow, aspow, nL, nH));
// }

// void Xt(double mb, double mW, double mZ, double mH, double mt, double mu, int apow, int aspow, int nL,int nH) 
// {
//   OSinput oi(mb, mW, mZ, mH, mt);

//   tt ttm = get_tt(oi, pow(mu,2));

//   MLPutReal128(stdlink, ttm.x(apow, aspow, nL, nH));
// }



// Pure QCD corrections

void XbQCD(double mb, double mW, double mZ, double mH, double mt, double mu, int nf) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  bb bbm = get_bb(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 3);

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xb", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 1);
  MLPutReal128(stdlink, bbm.x01(nf));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xb", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 2);
  MLPutReal128(stdlink, bbm.x02(nf));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xb", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 3);
  MLPutReal128(stdlink, bbm.x03(nf));

}

void XtQCD(double mb, double mW, double mZ, double mH, double mt, double mu, int nf) 
{
  OSinput oi(mb, mW, mZ, mH, mt);
  
  tt ttm = get_tt(oi, pow(mu,2));

  MLPutFunction(stdlink, "List", 3);

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xt", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 1);
  MLPutReal128(stdlink, ttm.x01(nf));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xt", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 2);
  MLPutReal128(stdlink, ttm.x02(nf));

  MLPutFunction(stdlink, "Rule", 2);
  MLPutFunction(stdlink, "xt", 2);
  MLPutInteger(stdlink, 0);
  MLPutInteger(stdlink, 3);
  MLPutReal128(stdlink, ttm.x03(nf));

}

int main(int argc, char* argv[]) 
{
  return MLMain(argc, argv);
}
