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

#ifndef __P2MS_HPP__
#define __P2MS_HPP__

#include "boost/math/tools/roots.hpp"
#include "boost/numeric/interval.hpp"

#include "tdecl.hpp"
#include "sminput.hpp"

#include "WW.hpp"
#include "ZZ.hpp"
#include "HH.hpp"
#include "tt.hpp"
#include "bb.hpp"
#include "dr.hpp"
#include "alphaGF.hpp"

namespace mr
{
  class AlphaBase
  {
  

  protected:
    OSinput oi;
    unsigned ord;
    double Gf0;
    double alphaS;
    
    double alpha0;

    Tolerance tol;
    static Tolerance getTolerance(double tol_)
    {
      return Tolerance(tol_);
    }

  public:
    AlphaBase(const OSinput & in_, double tol_ = 10e-9,
              const double &  Gf0_ = pdg2014::Gf,
              const double &  as_ = pdg2014::asMZ,
              unsigned order_ = 
              order::x01|order::x10|order::x02|
              order::x11|order::x20|order::x03): oi(in_), tol(getTolerance(tol_)), Gf0(Gf0_), alphaS(as_), ord(order_)
    {
     
    }

    virtual double operator()(const double& mu2 ) = 0;
    
  };


  // Alpha from solution of (eq.31) [hep-ph/1503.02138]
  class AlphaSolve : protected AlphaBase
  {
  
  public:
    AlphaSolve(const OSinput & in_, double tol_ = 10e-9,
               const double &  Gf0_ = pdg2014::Gf,
               const double &  as_ = pdg2014::asMZ,
               unsigned order_ = 
               order::x01|order::x10|order::x02|
               order::x11|order::x20|order::x03): AlphaBase(in_, tol_, Gf0_, as_,order_)
    {
     
    }

    double operator()(const double& mu2 );
    
  };

  // Alpha from solution of (eq.33) [hep-ph/1503.02138]
  class AlphaGF : protected AlphaBase
  {
  
public:
    AlphaGF(const OSinput & in_, double tol_ = 10e-9,
            const double &  Gf0_ = pdg2014::Gf,
            const double &  as_ = pdg2014::asMZ,
            unsigned order_ = 
            order::x01|order::x10|order::x02|
            order::x11|order::x20|order::x03): AlphaBase(in_, tol_, Gf0_, as_,order_)
    {
    
    }

    double operator()(const double& mu2 );
  };



  // SM version with fixed nL=2, nH=1
  template<class AlphaT>
  class P2MS
  {
    OSinput       oi;
    
    double  aEW;
    double aQCD;
    double   Gf;
    double   mu;

    // corrections
    double dbplus1;
    double dWplus1;
    double dZplus1;
    double dHplus1;
    double dtplus1;
    double dRplus1;

    unsigned ord;

  public:
    P2MS(const OSinput & oi_,
         const double &  Gf_ = pdg2014::Gf,
         const double &  as_ = pdg2014::asMZ,
         const double &  mu_ = pdg2014::MZ,
         unsigned ord_ = order::x01|order::x10|order::x02|order::x11|order::x20|order::x03 );
  
    double   a1() const;
    double   a2() const;
    double   as() const;
    double   at() const;
    double   ab() const;
    double alam() const;
  
  
    double  g1()  const;
    double  g2()  const;
    double  gs()  const;
    double  yt()  const;
    double  yb()  const;
    double lam()  const;
    double mphi() const;
    double vev()  const;
  
    MSinput getMSpar();

    std::vector<double> runningCouplings() const;

    SMCouplings ai() const;
    
    double scale() const
    {
      return mu;
    }
  
  };
  
  // SM NH and NL fixed
  template<class AlphaT>
  P2MS<AlphaT>::P2MS(const OSinput & oi_, const double &  Gf_, const double &  as_,const double &  mu_, unsigned ord_): oi(oi_), Gf(Gf_), aQCD(as_/4./Pi), mu(mu_), ord(ord_)
  {
    double mu2 = pow(mu,2);
  
    bb<OS> ibb(oi, mu2);
    WW<OS> iww(oi, mu2);
    ZZ<OS> izz(oi, mu2);
    HH<OS> ihh(oi, mu2);
    tt<OS> itt(oi, mu2);
    dr<OS> idr(oi, mu2);

    AlphaT alpha(oi, 10e-9, Gf, as_, ord);
  
    aEW = alpha(mu2)/4./Pi;
 
    const size_t fw = 4;
    lout(logINFO) << "\t  ---------------------------------------------------------------- ";
    lout(logINFO) << "\t  |                     Enabled corrections                      | ";
    lout(logINFO) << "\t  ---------------------------------------------------------------- ";
    lout(logINFO) << "\t  |   QCD  |   EW   |  QCD^2 | EW*QCD |  EW^2  |  QCD^3 |  QCD^4 | ";
    lout(logINFO) << "\t  |"
                  << std::setw(fw) <<  std::internal << (bool(ord & order::x01)? '+' : '-') << "    |"
                  << std::setw(fw) <<  std::internal << (bool(ord & order::x10)? '+' : '-') << "    |"
                  << std::setw(fw) <<  std::internal << (bool(ord & order::x02)? '+' : '-') << "    |"
                  << std::setw(fw) <<  std::internal << (bool(ord & order::x11)? '+' : '-') << "    |"
                  << std::setw(fw) <<  std::internal << (bool(ord & order::x20)? '+' : '-') << "    |"
                  << std::setw(fw) <<  std::internal << (bool(ord & order::x03)? '+' : '-') << "    |"
                  << std::setw(fw) <<  std::internal << (bool(ord & order::x04)? '+' : '-') << "    |";
    lout(logINFO) << "\t  ---------------------------------------------------------------- ";
  
    lout(logINFO) << "Using 1/alpha = " << 1./(4.*Pi*aEW);
    lout(logINFO) << "Matching scale mu = " << mu; 

    dbplus1 = 1;
    dWplus1 = 1;
    dZplus1 = 1;
    dHplus1 = 1;
    dtplus1 = 1;
    // And running vev
    dRplus1 = 1;



    if(ord & order::x01)
      {
        dbplus1 += aQCD*ibb.y01();
        dtplus1 += aQCD*itt.y01();
      }
    if(ord & order::x10)
      {
        dbplus1 += aEW*ibb.y10();
        dWplus1 += aEW*iww.y10();
        dZplus1 += aEW*izz.y10();
        dHplus1 += aEW*ihh.y10();
        dtplus1 += aEW*itt.y10();
        dRplus1 += aEW*idr.dr10();
      }
    if(ord & order::x02)
      {
        dbplus1 += aQCD*aQCD*ibb.y02();
        dtplus1 += aQCD*aQCD*itt.y02();
      }
    if(ord & order::x11)
      {
        dbplus1 += aEW*aQCD*ibb.y11();
        dWplus1 += aEW*aQCD*iww.y11();
        dZplus1 += aEW*aQCD*izz.y11();
        dHplus1 += aEW*aQCD*ihh.y11();
        dtplus1 += aEW*aQCD*itt.y11();
        dRplus1 += aEW*aQCD*idr.dr11();
      }
    if(ord & order::x20)
      {
        dbplus1 += aEW*aEW*ibb.y20();
        dWplus1 += aEW*aEW*iww.y20();
        dZplus1 += aEW*aEW*izz.y20();
        dHplus1 += aEW*aEW*ihh.y20();
        dtplus1 += aEW*aEW*itt.y20();
        dRplus1 += aEW*aEW*idr.dr20();

      }
    if(ord & order::x03)
      {
        dbplus1 += aQCD*aQCD*aQCD*ibb.y03();
        dtplus1 += aQCD*aQCD*aQCD*itt.y03();
      
      }
    if(ord & order::x04)
      {
        dbplus1 += aQCD*aQCD*aQCD*aQCD*ibb.y04();
        dtplus1 += aQCD*aQCD*aQCD*aQCD*itt.y04();
      
      }
  
    double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
    double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
  
    double a1 = 5./3.*(gg_ggp - gg)/16./Pi/Pi;
    double a2 = gg/16./Pi/Pi;
    double aS = aQCD;
    double ayt = pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
    double alam = Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;


    double vev2 = dRplus1/Gf/sqrt(2.);
  
    lout(logDEBUG) << " At matching scale mu = " << mu ;
    lout(logDEBUG) << " g1 = " << sqrt(3./5.*a1)*4*Pi;
    lout(logDEBUG) << " g2 = " << sqrt(a2)*4*Pi;
    lout(logDEBUG) << " g3 = " << sqrt(aS)*4*Pi;
    lout(logDEBUG) << " yt = " << sqrt(ayt)*4*Pi;
    lout(logDEBUG) << " lam = " << alam*16*Pi*Pi;
  
    lout(logDEBUG) << " vev = " << sqrt(vev2);
    
    lout(logDEBUG) << " mu0 = " << sqrt(2.*lam())*vev();

  }

  template<class AlphaT>
  double P2MS<AlphaT>::a1() const
  {
    double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
    double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
  
    return 5./3.*(gg_ggp - gg)/16./Pi/Pi;
  }

  template<class AlphaT>
  double P2MS<AlphaT>::g1() const
  {
    return sqrt(a1())*4*Pi;
  }


  template<class AlphaT>
  double P2MS<AlphaT>::a2() const
  {
    double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
    return gg/16./Pi/Pi;
  }

  template<class AlphaT>
  double P2MS<AlphaT>::g2() const
  {
    return sqrt(a2())*4*Pi;
  }

  template<class AlphaT>
  double P2MS<AlphaT>::as() const
  {
    return aQCD;
  }

  template<class AlphaT>
  double P2MS<AlphaT>::gs() const
  {
    return sqrt(as())*4*Pi;
  }
  
  template<class AlphaT>
  double P2MS<AlphaT>::at() const
  {
    return pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
  }

  template<class AlphaT>
  double P2MS<AlphaT>::yt() const
  {
    return sqrt(at())*4*Pi;  
  }


  template<class AlphaT>
  double P2MS<AlphaT>::ab() const
  {
    return pow(2.,3./2.)*Gf*oi.MMb()*pow(dbplus1,2)/16./Pi/Pi;
  }

  template<class AlphaT>
  double P2MS<AlphaT>::yb() const
  {
    return sqrt(ab())*4*Pi;  
  }


  template<class AlphaT>
  double P2MS<AlphaT>::alam() const
  {
    return Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;
  }

  template<class AlphaT>
  double P2MS<AlphaT>::lam() const
  {
    return alam()*16*Pi*Pi;
  }


  template<class AlphaT>
  double P2MS<AlphaT>::mphi() const // tree: mu0=Mh
  {
    return sqrt(2.*lam())*vev();
  }


  template<class AlphaT>
  double P2MS<AlphaT>::vev() const
  {
    return sqrt(dRplus1/Gf/sqrt(2.));
  }



  // MS input for conversion OS -> MS

  template<class AlphaT>
  MSinput P2MS<AlphaT>::getMSpar()
  {
  
    return MSinput::fromCouplings(g1(),     // U(1) 3/5
                                  g2(),     // SU(2) 
                                  yb(), 
                                  yt(), 
                                  lam(), 
                                  mphi(),   // Higgs mass parameter
                                           // normalized as mu0=Mh at
                                           // tree level
                                  mu);     // Input scale

  }

  template<class AlphaT>
  std::vector<double> P2MS<AlphaT>::runningCouplings() const
  {
    std::vector<double> g(9);

    g[couplings::g1] = g1();
    g[couplings::g2] = g2();
    g[couplings::gs] = gs();
    g[couplings::yt] = yt();
    g[couplings::yb] = yb();
    g[couplings::ytau] = 0;
    g[couplings::lam] = lam();
    g[couplings::mphi] = mphi();
    g[couplings::vev] = vev();

    return g;
  }


  template<class AlphaT>
  SMCouplings P2MS<AlphaT>::ai() const
  {
    SMCouplings a(9);

    a[couplings::g1] = a1();
    a[couplings::g2] = a2();
    a[couplings::gs] = as();
    a[couplings::yt] = at();
    a[couplings::yb] = ab();
    a[couplings::ytau] = 0;
    a[couplings::lam] = alam();
    a[couplings::mphi] = mphi();
    a[couplings::vev] = vev();

    return a;
  }


} // namespace mr


#endif  // __P2MS_HPP__

