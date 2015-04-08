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

#include "boost/math/tools/roots.hpp"
#include "boost/numeric/interval.hpp"
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
  unsigned ord;
  long double Gf0;
  long double alphaS;
  
  DiffGF(OSinput in_, long double Gf0_, long double as_, long double mu2, unsigned order_) : oi(in_), Gf0(Gf0_), alphaS(as_), ord(order_)
  {
    dW = new WW<OS>(oi, mu2);
    dZ = new ZZ<OS>(oi, mu2);
  }
  
  long double operator()(long double alpha)
  {
    long double dMyW = 1;
    long double dMyZ = 1;
    long double Gf;
    
    if(ord & order::x10)
      {
        dMyW += alpha/4./Pi*dW->y10();
        dMyZ += alpha/4./Pi*dZ->y10();
      }
    if(ord & order::x11)
      {
        dMyW += alpha/4./Pi*alphaS/4./Pi*dW->y11();
        dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ->y11();
      }
    if(ord & order::x20)
      {
        dMyW += pow(alpha/4./Pi,2)*dW->y20();
        dMyZ += pow(alpha/4./Pi,2)*dZ->y20();
      }
    

    Gf = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
    return (Gf0 - Gf)*pow(10,2);
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




P2MSnLnH::P2MSnLnH(const OSinput & oi_, const long double &  Gf_, const long double &  as_,const long double &  mu_): oi(oi_), Gf(Gf_), aQCD(as_/4./Pi), mu(mu_)
{
  long double mu2 = pow(mu,2);

  bp  = new bb<OS>(oi, mu2);
  wp  = new WW<OS>(oi, mu2);
  zp  = new ZZ<OS>(oi, mu2);
  hp  = new HH<OS>(oi, mu2);
  tp  = new tt<OS>(oi, mu2);
  drp = new dr<OS>(oi, mu2);


  DiffGF dGF(oi, Gf, as_, mu2, order::all);
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


long double P2MSnLnH::a1(size_t nL, size_t nH)
{
  long double dWplus1 = 1 + aEW*wp->y10(nL, nH) + aEW*aQCD*wp->y11(nL, nH) + aEW*aEW*wp->y20(nL, nH);
  long double dZplus1 = 1 + aEW*zp->y10(nL, nH) + aEW*aQCD*zp->y11(nL, nH) + aEW*aEW*zp->y20(nL, nH);
  
  long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
  long double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
  
  return 5./3.*(gg_ggp - gg)/16./Pi/Pi;
  
  
}

long double P2MSnLnH::g1(size_t nL, size_t nH)
{
  return sqrt(3./5.*a1(nL,nH))*4*Pi;
}



long double P2MSnLnH::a2(size_t nL, size_t nH)
{
  long double dWplus1 = 1 + aEW*wp->y10(nL, nH) + aEW*aQCD*wp->y11(nL, nH) + aEW*aEW*wp->y20(nL, nH);

  long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
  
  return gg/16./Pi/Pi;
}

long double P2MSnLnH::g2(size_t nL, size_t nH)
{
  return sqrt(a2(nL,nH))*4*Pi;
}



long double P2MSnLnH::as(size_t nL, size_t nH)
{
  return aQCD;
}

long double P2MSnLnH::gs(size_t nL, size_t nH)
{
  return sqrt(as())*4*Pi;
}



long double P2MSnLnH::at(size_t nL, size_t nH)
{
  long double dtplus1 = 1 + aEW*tp->y10(nL, nH) + aEW*aQCD*tp->y11(nL, nH) + aEW*aEW*tp->y20(nL, nH)
    + aQCD*tp->y01(nL, nH)+ aQCD*aQCD*tp->y02(nL, nH)+ aQCD*aQCD*aQCD*tp->y03(nL, nH);
  
  return pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
}

long double P2MSnLnH::yt(size_t nL, size_t nH)
{
  return sqrt(at(nL,nH))*4*Pi;  
}



long double P2MSnLnH::ab(size_t nL, size_t nH)
{
  long double dbplus1 = 1 + aEW*bp->y10(nL, nH) + aEW*aQCD*bp->y11(nL, nH) + aEW*aEW*bp->y20(nL, nH)
    + aQCD*bp->y01(nL, nH)+ aQCD*aQCD*bp->y02(nL, nH)+ aQCD*aQCD*aQCD*bp->y03(nL, nH);
                                                                    
  return pow(2.,3./2.)*Gf*oi.MMb()*pow(dbplus1,2)/16./Pi/Pi;
}

long double P2MSnLnH::yb(size_t nL, size_t nH)
{
  return sqrt(ab(nL,nH))*4*Pi;  
}



long double P2MSnLnH::alam(size_t nL, size_t nH)
{
  long double dHplus1 = 1 + aEW*hp->y10(nL, nH) + aEW*aQCD*hp->y11(nL, nH) + aEW*aEW*hp->y20(nL, nH);
  
  return Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;
}

long double P2MSnLnH::lam(size_t nL, size_t nH)
{
  return alam(nL,nH)*16*Pi*Pi;
}



long double P2MSnLnH::mu0(size_t nL, size_t nH) // tree: mu0=Mh
{
  return sqrt(2.*lam())*vev();
}



long double P2MSnLnH::vev(size_t nL, size_t nH)
{
  long double dRplus1 = 1 + aEW*drp->dr10() + aEW*aQCD*drp->dr11() + aEW*aEW*drp->dr20();

  return sqrt(dRplus1/Gf/sqrt(2.));
}



// MS input for conversion OS -> MS

MSinput P2MSnLnH::getMSpar()
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

std::vector<long double> P2MSnLnH::runningCouplings()
{
  std::vector<long double> g(9);

  g[couplings::g1] = g1();
  g[couplings::g2] = g2();
  g[couplings::gs] = gs();
  g[couplings::yt] = yt();
  g[couplings::yb] = yb();
  g[couplings::ytau] = 0;
  g[couplings::lam] = lam();
  g[couplings::mu0] = mu0();
  g[couplings::vev] = vev();

  return g;
}



// Simplified version
P2MS::P2MS(const OSinput & oi_, const long double &  Gf_, const long double &  as_,const long double &  mu_, unsigned ord_): oi(oi_), Gf(Gf_), aQCD(as_/4./Pi), mu(mu_), ord(ord_)
{
  long double mu2 = pow(mu,2);
  
  bp  = new bb<OS>(oi, mu2);
  wp  = new WW<OS>(oi, mu2);
  zp  = new ZZ<OS>(oi, mu2);
  hp  = new HH<OS>(oi, mu2);
  tp  = new tt<OS>(oi, mu2);
  drp = new dr<OS>(oi, mu2);
  

  DiffGF dGF(oi, Gf, as_, mu2, ord);
  tolerance tol = 1e-12;
  
  std::pair<long double, long double> found = boost::math::tools::bisect(dGF, 1./140., 1./120., tol);
  
  boost::numeric::interval<long double> fint(found.first, found.second);
  
  std::cout << std::setprecision(10);
  std::cout << "==> 1/alpha = [" << 1./found.first << ',' << 1./found.second << "]\n";
  

  aEW  = boost::numeric::median(fint)/4./Pi;

  std::cout << "alpha/4/pi = " << aEW << std::endl;


  std::cout << "Enabled corrections:" << std::endl;

  const size_t fw = 4;
  std::cout << "\t  |   QCD  |   EW   |  QCD^2 | EW*QCD |  EW^2  |  QCD^3 | " << std::endl;
  std::cout << "\t  |"
            << std::setw(fw) <<  std::internal << (bool(ord & order::x01)? '+' : '-') << "    |"
            << std::setw(fw) <<  std::internal << bool(ord & order::x10) << "    |"
            << std::setw(fw) <<  std::internal << bool(ord & order::x02) << "    |"
            << std::setw(fw) <<  std::internal << bool(ord & order::x11) << "    |"
            << std::setw(fw) <<  std::internal << bool(ord & order::x20) << "    |"
            << std::setw(fw) <<  std::internal << bool(ord & order::x03) << "    |"
            << std::endl;

  
  

  dbplus1 = 1;
  dWplus1 = 1;
  dZplus1 = 1;
  dHplus1 = 1;
  dtplus1 = 1;
  // And running vev
  dRplus1 = 1;



  if(ord & order::x01)
    {
      dbplus1 += aQCD*bp->y01();
      dtplus1 += aQCD*tp->y01();
    }
  if(ord & order::x10)
    {
      dbplus1 += aEW*bp->y10();
      dWplus1 += aEW*wp->y10();
      dZplus1 += aEW*zp->y10();
      dHplus1 += aEW*hp->y10();
      dtplus1 += aEW*tp->y10();
      dRplus1 += aEW*drp->dr10();
    }
  if(ord & order::x02)
    {
      dbplus1 += aQCD*aQCD*bp->y02();
      dtplus1 += aQCD*aQCD*tp->y02();
    }
  if(ord & order::x11)
    {
      dbplus1 += aEW*aQCD*bp->y11();
      dWplus1 += aEW*aQCD*wp->y11();
      dZplus1 += aEW*aQCD*zp->y11();
      dHplus1 += aEW*aQCD*hp->y11();
      dtplus1 += aEW*aQCD*tp->y11();
      dRplus1 += aEW*aQCD*drp->dr11();
    }
  if(ord & order::x20)
    {
      dbplus1 += aEW*aEW*bp->y20();
      dWplus1 += aEW*aEW*wp->y20();
      dZplus1 += aEW*aEW*zp->y20();
      dHplus1 += aEW*aEW*hp->y20();
      dtplus1 += aEW*aEW*tp->y20();
      dRplus1 += aEW*aEW*drp->dr20();

    }
  if(ord & order::x03)
    {
      dbplus1 += aQCD*aQCD*aQCD*bp->y03();
      dtplus1 += aQCD*aQCD*aQCD*tp->y03();
      
    }
  

  
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


long double P2MS::a1() const
{
  long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
  long double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
  
  return 5./3.*(gg_ggp - gg)/16./Pi/Pi;
}

long double P2MS::g1() const
{
  return sqrt(3./5.*a1())*4*Pi;
}



long double P2MS::a2() const
{
  long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
  return gg/16./Pi/Pi;
}

long double P2MS::g2() const
{
  return sqrt(a2())*4*Pi;
}

long double P2MS::as() const
{
  return aQCD;
}

long double P2MS::gs() const
{
  return sqrt(as())*4*Pi;
}
  

long double P2MS::at() const
{
  return pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
}

long double P2MS::yt() const
{
  return sqrt(at())*4*Pi;  
}



long double P2MS::ab() const
{
  return pow(2.,3./2.)*Gf*oi.MMb()*pow(dbplus1,2)/16./Pi/Pi;
}

long double P2MS::yb() const
{
  return sqrt(ab())*4*Pi;  
}



long double P2MS::alam() const
{
  return Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;
}

long double P2MS::lam() const
{
  return alam()*16*Pi*Pi;
}



long double P2MS::mu0() const // tree: mu0=Mh
{
  return sqrt(2.*lam())*vev();
}



long double P2MS::vev() const
{
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

std::vector<long double> P2MS::runningCouplings() const
{
  std::vector<long double> g(9);

  g[couplings::g1] = g1();
  g[couplings::g2] = g2();
  g[couplings::gs] = gs();
  g[couplings::yt] = yt();
  g[couplings::yb] = yb();
  g[couplings::ytau] = 0;
  g[couplings::lam] = lam();
  g[couplings::mu0] = mu0();
  g[couplings::vev] = vev();

  return g;
}


std::vector<long double> P2MS::ai() const
{
  std::vector<long double> a(9);

  a[couplings::g1] = a1();
  a[couplings::g2] = a2();
  a[couplings::gs] = as();
  a[couplings::yt] = at();
  a[couplings::yb] = ab();
  a[couplings::ytau] = 0;
  a[couplings::lam] = alam();
  a[couplings::mu0] = mu0();
  a[couplings::vev] = vev();

  return a;
}
