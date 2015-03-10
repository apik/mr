//
// MR - 2-loop matching and 3-loop Running, including full 2-loop EW corrections
// Copyright (C) 2014 Andrey Pikelner <pikelner@theor.jinr.ru>
//
// This file is part of MR.
//
// MR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MR.  If not, see <http://www.gnu.org/licenses/>.
//

#include <boost/math/tools/roots.hpp>
#include <boost/numeric/interval.hpp>
#include "p2ms.hpp"
#include <stdexcept>
#include "betaQCD.hpp"
#include "betaQEDQCD.hpp"
#include "bb.hpp" 
#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp" 



struct DiffGF
{
  
  OSinput oi;
  WW<OS>* dW;
  ZZ<OS>* dZ;
  size_t order;
  long double Gf0;
  long double alphaS;
  
  DiffGF(OSinput in_, long double Gf0_, long double as_, long double mu2, size_t order_) : oi(in_), Gf0(Gf0_), alphaS(as_), order(order_)
  {
    dW = new WW<OS>(oi, mu2);
    dZ = new ZZ<OS>(oi, mu2);
  }
  
  long double operator()(long double alpha)
  {
    long double dMyW = 1;
    long double dMyZ = 1;
    long double Gf[4];
    
    // Tree level
    Gf[0] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);

    // 1-loop level
    dMyW += alpha/4./Pi*dW->y10();
    dMyZ += alpha/4./Pi*dZ->y10();
    
    Gf[1] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
    // 2-loop level
    dMyW += alpha/4./Pi*alphaS/4./Pi*dW->y11();
    dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ->y11();
    
    Gf[2] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
    
    dMyW += pow(alpha/4./Pi,2)*dW->y20();
    dMyZ += pow(alpha/4./Pi,2)*dZ->y20();
    
    Gf[3] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
    return (Gf0 - Gf[order])*pow(10,2);
  }

};


class tolerance {
public:
  tolerance(long double eps) :
    _eps(eps) {
  }
  bool operator()(long double a, long double b) {
    return (fabs(b - a) <= _eps);
  }
private:
  long double _eps;
};




P2MS::P2MS(const OSinput & oi_, const long double &  Gf_, const long double &  as_,const long double &  mu_): oi(oi_), Gf(Gf_), aQCD(as_/4./Pi), mu(mu_)
{
  long double mu2 = pow(mu,2);

  bp  = new bb(oi, mu2);
  wp  = new WW<OS>(oi, mu2);
  zp  = new ZZ<OS>(oi, mu2);
  hp  = new HH<OS>(oi, mu2);
  tp  = new tt<OS>(oi, mu2);
  drp = new dr<OS>(oi, mu2);


  DiffGF dGF(oi, Gf, as_, mu2, 3);
  tolerance tol = 1e-12;
  
  std::pair<long double, long double> found = boost::math::tools::bisect(dGF, 1./140., 1./120., tol);
  
  boost::numeric::interval<long double> fint(found.first, found.second);
  
  std::cout << std::setprecision(10);
  std::cout << "==> 1/alpha = [" << 1./found.first << ',' << 1./found.second << "]\n";
  

  aEW  = boost::numeric::median(fint)/4./Pi;

  std::cout << "alpha/4/pi = " << aEW << std::endl;

  long double dbplus1 = 1 + aEW*bp->y10() + aEW*aQCD*bp->y11() + aEW*aEW*bp->y20()
    + aQCD*bp->y01()+ aQCD*aQCD*bp->y02()+ aQCD*aQCD*aQCD*bp->y03();
  long double dWplus1 = 1 + aEW*wp->y10() + aEW*aQCD*wp->y11() + aEW*aEW*wp->y20();
  long double dZplus1 = 1 + aEW*zp->y10() + aEW*aQCD*zp->y11() + aEW*aEW*zp->y20();
  long double dHplus1 = 1 + aEW*hp->y10() + aEW*aQCD*hp->y11() + aEW*aEW*hp->y20();
  long double dtplus1 = 1 + aEW*tp->y10() + aEW*aQCD*tp->y11() + aEW*aEW*tp->y20()
    + aQCD*tp->y01()+ aQCD*aQCD*tp->y02()+ aQCD*aQCD*aQCD*tp->y03();
  // And running vev
  long double dRplus1 = 1 + aEW*drp->dr10() + aEW*aQCD*drp->dr11() + aEW*aEW*drp->dr20();

  long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
  long double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
  
  long double a1 = 5./3.*(gg_ggp - gg)/16./Pi/Pi;
  long double a2 = gg/16./Pi/Pi;
  long double aS = aQCD;
  long double ayt = pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
  long double alam = Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;


  long double vev2 = dRplus1/Gf/sqrt(2.);
  
  std::cout << " At matching scale mu = " << mu << std::endl;
  std::cout << " g1 = " << sqrt(3./5.*a1)*4*Pi << std::endl;
  std::cout << " g2 = " << sqrt(a2)*4*Pi << std::endl;
  std::cout << " g3 = " << sqrt(aS)*4*Pi << std::endl;
  std::cout << " yt = " << sqrt(ayt)*4*Pi << std::endl;
  std::cout << " lam = " << alam*16*Pi*Pi << std::endl;
  
  std::cout << " vev = " << sqrt(vev2) << std::endl;

  std::cout << " mu0 = " << sqrt(2.*lam())*vev() << std::endl;
   
}


long double P2MS::a1(size_t nL, size_t nH)
{
  long double dWplus1 = 1 + aEW*wp->y10(nL, nH) + aEW*aQCD*wp->y11(nL, nH) + aEW*aEW*wp->y20(nL, nH);
  long double dZplus1 = 1 + aEW*zp->y10(nL, nH) + aEW*aQCD*zp->y11(nL, nH) + aEW*aEW*zp->y20(nL, nH);
  
  long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
  long double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
  
  return 5./3.*(gg_ggp - gg)/16./Pi/Pi;
  
  
}

long double P2MS::g1(size_t nL, size_t nH)
{
  return sqrt(3./5.*a1(nL,nH))*4*Pi;
}



long double P2MS::a2(size_t nL, size_t nH)
{
  long double dWplus1 = 1 + aEW*wp->y10(nL, nH) + aEW*aQCD*wp->y11(nL, nH) + aEW*aEW*wp->y20(nL, nH);

  long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
  
  return gg/16./Pi/Pi;
}

long double P2MS::g2(size_t nL, size_t nH)
{
  return sqrt(a2(nL,nH))*4*Pi;
}
  

long double P2MS::at(size_t nL, size_t nH)
{
  long double dtplus1 = 1 + aEW*tp->y10(nL, nH) + aEW*aQCD*tp->y11(nL, nH) + aEW*aEW*tp->y20(nL, nH)
    + aQCD*tp->y01(nL, nH)+ aQCD*aQCD*tp->y02(nL, nH)+ aQCD*aQCD*aQCD*tp->y03(nL, nH);
  
  return pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
}

long double P2MS::yt(size_t nL, size_t nH)
{
  return sqrt(at(nL,nH))*4*Pi;  
}



long double P2MS::ab(size_t nL, size_t nH)
{
  long double dbplus1 = 1 + aEW*bp->y10(nL, nH) + aEW*aQCD*bp->y11(nL, nH) + aEW*aEW*bp->y20(nL, nH)
    + aQCD*bp->y01(nL, nH)+ aQCD*aQCD*bp->y02(nL, nH)+ aQCD*aQCD*aQCD*bp->y03(nL, nH);
                                                                    
  return pow(2.,3./2.)*Gf*oi.MMt()*pow(dbplus1,2)/16./Pi/Pi;
}

long double P2MS::yb(size_t nL, size_t nH)
{
  return sqrt(ab(nL,nH))*4*Pi;  
}



long double P2MS::alam(size_t nL, size_t nH)
{
  long double dHplus1 = 1 + aEW*hp->y10(nL, nH) + aEW*aQCD*hp->y11(nL, nH) + aEW*aEW*hp->y20(nL, nH);
  
  return Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;
}

long double P2MS::lam(size_t nL, size_t nH)
{
  return alam(nL,nH)*16*Pi*Pi;
}



long double P2MS::mu0(size_t nL, size_t nH) // tree: mu0=Mh
{
  return sqrt(2.*lam())*vev();
}



long double P2MS::vev(size_t nL, size_t nH)
{
  long double dRplus1 = 1 + aEW*drp->dr10() + aEW*aQCD*drp->dr11() + aEW*aEW*drp->dr20();

  return sqrt(dRplus1/Gf/sqrt(2.));
}



// MS input for conversion OS -> MS

MSinput P2MS::getMSpar()
{
  
  return MSinput::fromConsts(mu, // Input scale
                             mu0(),   //Higgs mass parameter
                             //normalized as mu0=Mh at
                             //tree level
                             lam(), 
                             yb(), 
                             yt(), 
                             g2(),     // SU(2) 
                             g1()     // U(1)
                             );
}
