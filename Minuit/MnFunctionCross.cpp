#include "Minuit/MnFunctionCross.h"
#include "Minuit/FunctionMinimum.h"
#include "Minuit/MnMigrad.h"
#include "Minuit/FCNBase.h"
#include "Minuit/MnParabola.h"
#include "Minuit/MnParabolaPoint.h"
#include "Minuit/MnParabolaFactory.h"
#include "Minuit/MnCross.h"
#include "Minuit/MnMachinePrecision.h"
#include "Minuit/MnPrint.h"

#define DEBUG

MnCross MnFunctionCross::operator()(const std::vector<unsigned int>& par, const std::vector<double>& pmid, const std::vector<double>& pdir, double tlr, unsigned int maxcalls) const {
//   double edmmax = 0.5*0.001*toler*theFCN.up();
  
  //std::cout<<"MnFunctionCross par "<<par.front()<<std::endl;

  unsigned int npar = par.size();
  unsigned int nfcn = 0;
  const MnMachinePrecision& prec = theState.precision();
//   double tlr = 0.01;
  double tlf = tlr*theFCN.up();
  double tla = tlr;
  unsigned int maxitr = 15;
  unsigned int ipt = 0;
  double aminsv = theFval;
  double aim = aminsv + theFCN.up();
  //std::cout<<"aim= "<<aim<<std::endl;
  double aopt = 0.;
  bool limset = false;
  std::vector<double> alsb(3, 0.), flsb(3, 0.);
  double up = theFCN.up();

  double aulim = 100.;
  for(unsigned int i = 0; i < par.size(); i++) {
    unsigned int kex = par[i];
    if(theState.parameter(kex).hasLimits()) {
      double zmid = pmid[i];
      double zdir = pdir[i];
      if(fabs(zdir) < theState.precision().eps()) continue;
//       double zlim = 0.;
      if(zdir > 0. && theState.parameter(kex).hasUpperLimit()) {
	double zlim = theState.parameter(kex).upperLimit();
	aulim = std::min(aulim, (zlim-zmid)/zdir);
      }
      else if(zdir < 0. && theState.parameter(kex).hasLowerLimit()) {
	double zlim = theState.parameter(kex).lowerLimit();
	aulim = std::min(aulim, (zlim-zmid)/zdir);
      }
    }
  }
  //std::cout<<"100"<<std::endl;

  if(aulim  < aopt+tla) limset = true;

  MnMigrad migrad(theFCN, theState, MnStrategy(std::max(0, int(theStrategy.strategy()-1))));

  for(unsigned int i = 0; i < npar; i++) {
    migrad.setValue(par[i], pmid[i]);
  }

  FunctionMinimum min0 = migrad(maxcalls, tlr);
  nfcn += min0.nfcn();

  if(min0.hasReachedCallLimit()) 
    return MnCross(min0.userState(), nfcn, MnCross::CrossFcnLimit());
  if(!min0.isValid()) return MnCross(theState, nfcn);
  if(limset == true && min0.fval() < aim) 
    return MnCross(min0.userState(), nfcn, MnCross::CrossParLimit());

  ipt++;
  alsb[0] = 0.;
  flsb[0] = min0.fval();
  flsb[0] = std::max(flsb[0], aminsv + 0.1*up);
  aopt = sqrt(up/(flsb[0]-aminsv)) - 1.;
  if(fabs(flsb[0] - aim) < tlf) return MnCross(aopt, min0.userState(), nfcn);

  if(aopt > 1.) aopt = 1.;
  if(aopt < -0.5) aopt = -0.5;
  limset = false;
  if(aopt > aulim) {
    aopt = aulim;
    limset = true;
  }

  for(unsigned int i = 0; i < npar; i++) {
    migrad.setValue(par[i], pmid[i] + (aopt)*pdir[i]);
  }

  FunctionMinimum min1 = migrad(maxcalls, tlr);
  nfcn += min1.nfcn();

  if(min1.hasReachedCallLimit()) 
    return MnCross(min1.userState(), nfcn, MnCross::CrossFcnLimit());
  if(!min1.isValid()) return MnCross(theState, nfcn);
  if(limset == true && min1.fval() < aim) 
    return MnCross(min1.userState(), nfcn, MnCross::CrossParLimit());
  
  ipt++;
  alsb[1] = aopt;
  flsb[1] = min1.fval();
  double dfda = (flsb[1] - flsb[0])/(alsb[1] - alsb[0]);
//   if(dfda > 0.) goto L460;

L300:
  if(dfda < 0.) {
    //std::cout<<"300"<<std::endl;
    unsigned int maxlk = maxitr - ipt;
    for(unsigned int it = 0; it < maxlk; it++) {
      alsb[0] = alsb[1];
      flsb[0] = flsb[1];
      // LM: add + 1, looking at Fortran code it starts from 1 ( see bug #8396)
      aopt = alsb[0] + 0.2*(it+1);   
      limset = false;
      if(aopt > aulim) {
	aopt = aulim;
	limset = true;
      }
      for(unsigned int i = 0; i < npar; i++) {
	migrad.setValue(par[i], pmid[i] + (aopt)*pdir[i]);
      }
      min1 = migrad(maxcalls, tlr);
      nfcn += min1.nfcn();

      if(min1.hasReachedCallLimit()) 
	return MnCross(min1.userState(), nfcn, MnCross::CrossFcnLimit());
      if(!min1.isValid()) return MnCross(theState, nfcn);
      if(limset == true && min1.fval() < aim) 
	return MnCross(min1.userState(), nfcn, MnCross::CrossParLimit());
      ipt++;
      alsb[1] = aopt;
      flsb[1] = min1.fval();
      dfda = (flsb[1] - flsb[0])/(alsb[1] - alsb[0]);
//       if(dfda > 0.) goto L460;
      if(dfda > 0.) break;
    } 
    if(ipt > maxitr) return MnCross(theState, nfcn);
  } //if(dfda < 0.)
  
L460:
  //std::cout<<"460"<<std::endl;
    
  aopt = alsb[1] + (aim-flsb[1])/dfda;
  double fdist = std::min(fabs(aim  - flsb[0]), fabs(aim  - flsb[1]));
  double adist = std::min(fabs(aopt - alsb[0]), fabs(aopt - alsb[1]));
  tla = tlr;
  if(fabs(aopt) > 1.) tla = tlr*fabs(aopt);
  if(adist < tla && fdist < tlf) return MnCross(aopt, min1.userState(), nfcn);
  if(ipt > maxitr) return MnCross(theState, nfcn);
  double bmin = std::min(alsb[0], alsb[1]) - 1.;
  if(aopt < bmin) aopt = bmin;
  double bmax = std::max(alsb[0], alsb[1]) + 1.;
  if(aopt > bmax) aopt = bmax;
  
  limset = false;
  if(aopt > aulim) {
    aopt = aulim;
    limset = true;
  }

  for(unsigned int i = 0; i < npar; i++) {
    migrad.setValue(par[i], pmid[i] + (aopt)*pdir[i]);
  }
  FunctionMinimum min2 = migrad(maxcalls, tlr);
  nfcn += min2.nfcn();

  if(min2.hasReachedCallLimit()) 
    return MnCross(min2.userState(), nfcn, MnCross::CrossFcnLimit());
  if(!min2.isValid()) return MnCross(theState, nfcn);
  if(limset == true && min2.fval() < aim) 
    return MnCross(min2.userState(), nfcn, MnCross::CrossParLimit());

  ipt++;
  alsb[2] = aopt;
  flsb[2] = min2.fval();

  double ecarmn = fabs(flsb[2] - aim);
  double ecarmx = 0.;
  unsigned int ibest = 2;
  unsigned int iworst = 0;
  unsigned int noless = 0;

  for(unsigned int i = 0; i < 3; i++) {
    double ecart = fabs(flsb[i] - aim);
    if(ecart > ecarmx) {
      ecarmx = ecart;
      iworst = i;
    }
    if(ecart < ecarmn) {
      ecarmn = ecart;
      ibest = i;
    }
    if(flsb[i] < aim) noless++;
  }

  //std::cout<<"480"<<std::endl;

  if(noless == 1 || noless == 2) goto L500;
  if(noless == 0 && ibest != 2) return MnCross(theState, nfcn);
  if(noless == 3 && ibest != 2) {
    alsb[1] = alsb[2];
    flsb[1] = flsb[2];
    goto L300;
  }

  flsb[iworst] = flsb[2];
  alsb[iworst] = alsb[2];
  dfda = (flsb[1] - flsb[0])/(alsb[1] - alsb[0]);
  goto L460;

L500:
  //std::cout<<"500"<<std::endl;
  do {
    MnParabola parbol = MnParabolaFactory()(MnParabolaPoint(alsb[0], flsb[0]), MnParabolaPoint(alsb[1], flsb[1]), MnParabolaPoint(alsb[2], flsb[2]));
    //   aopt = parbol.x_pos(aim);
    //std::cout<<"alsb1,2,3= "<<alsb[0]<<", "<<alsb[1]<<", "<<alsb[2]<<std::endl;
    //std::cout<<"flsb1,2,3= "<<flsb[0]<<", "<<flsb[1]<<", "<<flsb[2]<<std::endl;
    
    double coeff1 = parbol.c();
    double coeff2 = parbol.b(); 
    double coeff3 = parbol.a(); 
    double determ = coeff2*coeff2 - 4.*coeff3*(coeff1 - aim);
    if(determ < prec.eps()) return MnCross(theState, nfcn); 
    double rt = sqrt(determ);
    double x1 = (-coeff2 + rt)/(2.*coeff3);
    double x2 = (-coeff2 - rt)/(2.*coeff3);
    double s1 = coeff2 + 2.*x1*coeff3;
    double s2 = coeff2 + 2.*x2*coeff3;
    
    if(s1*s2 > 0.) std::cout<<"MnFunctionCross problem 1"<<std::endl;
    aopt = x1;
    double slope = s1;
    if(s2 > 0.) {
      aopt = x2;
      slope = s2;
    }
    
    tla = tlr;
    if(fabs(aopt) > 1.) tla = tlr*fabs(aopt);
    if(fabs(aopt - alsb[ibest]) < tla && fabs(flsb[ibest] - aim) < tlf) 
      return MnCross(aopt, min2.userState(), nfcn);
    
//     if(ipt > maxitr) return MnCross();
    
    unsigned int ileft = 3;
    unsigned int iright = 3;
    unsigned int iout = 3;
    ibest = 0;
    ecarmx = 0.;
    ecarmn = fabs(aim-flsb[0]);
    for(unsigned int i = 0; i < 3; i++) {
      double ecart = fabs(flsb[i] - aim);
      if(ecart < ecarmn) {
	ecarmn = ecart;
	ibest = i;
      }
      if(ecart > ecarmx) ecarmx = ecart;
      if(flsb[i] > aim) {
	if(iright == 3) iright = i;
	else if(flsb[i] > flsb[iright]) iout = i;
	else {
	  iout = iright;
	  iright = i;
	}
      } else if(ileft == 3) ileft = i;
      else if(flsb[i] < flsb[ileft]) iout = i;
      else {
	iout = ileft;
	ileft = i;
      }
    }
    //std::cout<<"550"<<std::endl;
    
    if(ecarmx > 10.*fabs(flsb[iout] - aim)) 
      aopt = 0.5*(aopt + 0.5*(alsb[iright] + alsb[ileft]));
    double smalla = 0.1*tla;
    if(slope*smalla > tlf) smalla = tlf/slope;
    double aleft = alsb[ileft] + smalla;
    double aright = alsb[iright] - smalla;
    if(aopt < aleft) aopt = aleft;
    if(aopt > aright) aopt = aright;
    if(aleft > aright) aopt = 0.5*(aleft + aright);
    
    limset = false;
    if(aopt > aulim) {
      aopt = aulim;
      limset = true;
    }
    
    for(unsigned int i = 0; i < npar; i++) {
      migrad.setValue(par[i], pmid[i] + (aopt)*pdir[i]);
    }
    min2 = migrad(maxcalls, tlr);
    nfcn += min2.nfcn();

    if(min2.hasReachedCallLimit()) 
      return MnCross(min2.userState(), nfcn, MnCross::CrossFcnLimit());
    if(!min2.isValid()) return MnCross(nfcn);
    if(limset == true && min2.fval() < aim) 
      return MnCross(min2.userState(), nfcn, MnCross::CrossParLimit());
    
    ipt++;
    alsb[iout] = aopt;
    flsb[iout] = min2.fval();
    ibest = iout;
  } while(ipt < maxitr);
  
  // goto L500;

  return MnCross(theState, nfcn);
}
