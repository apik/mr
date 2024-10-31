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

#ifndef __SMINPUT_HPP__
#define __SMINPUT_HPP__
#include "constants.hpp"
#include <stdexcept>

namespace mr
{
  template<class T>
  class OSinputTemplate
  {
    T iMb;
    T iMW;
    T iMZ;
    T iMH;
    T iMt;

  public:
    OSinputTemplate(T Mb_, T MW_, T MZ_, T MH_, T Mt_) : 
      iMb(Mb_), iMW(MW_), iMZ(MZ_), iMH(MH_), iMt(Mt_)
    {
    }

    OSinputTemplate() 
    {
    }
  
    bool operator < (const OSinputTemplate& b) const 
    {
      if(iMb < b.Mb()) return true;
      else if(iMb > b.Mb()) return false;
      else
      
        if(iMW < b.MW()) return true;
        else if(iMW > b.MW()) return false;
        else
        
          if(iMZ < b.MZ()) return true;
          else if(iMZ > b.MZ()) return false;
          else
        
            if(iMH < b.MH()) return true;
            else if(iMH > b.MH()) return false;
            else
              if(iMt < b.Mt()) return true;
              else if(iMt > b.Mt()) return false;
      return false;
    }

    T MMb() const
    {
      return iMb*iMb;
    }
    T MMW() const
    {
      return iMW*iMW;
    }
    T MMZ() const
    {
      return iMZ*iMZ;
    }
    T MMH() const
    {
      return iMH*iMH;
    }
    T MMt() const
    {
      return iMt*iMt;
    }

    // Weinberg trigonometric
    T CW() const
    {
      return iMW/iMZ;
    }
    T CCW() const
    {
      return pow(iMW/iMZ,2);
    }
    T SW() const
    {
      return sqrt(1-CCW());
    }
    T SSW() const
    {
      return 1-CCW();
    }


    T Mb() const
    {
      return iMb;
    }
    T MW() const
    {
      return iMW;
    }
    T MZ() const
    {
      return iMZ;
    }
    T MH() const
    {
      return iMH;
    }
    T Mt() const
    {
      return iMt;
    }

    // Modification
    OSinputTemplate<T>& setMb(T mb)
    {
      iMb = mb;
      return *this;
    }
    OSinputTemplate<T>& setMW(T mW)
    {
      iMW = mW;
      return *this;
    }
    OSinputTemplate<T>& setMZ(T mZ)
    {
      iMZ = mZ;
      return *this;
    }
    OSinputTemplate<T>& setMH(T mH)
    {
      iMH = mH;
      return *this;
    }
    OSinputTemplate<T>& setMt(T mt)
    {
      iMt = mt;
      return *this;
    }

  };

  typedef OSinputTemplate<double> OSinput;



  class MSinput
  {
    double imb;
    double imW;
    double imZ;
    double imH;
    double imt;
    // only MS
    double iv;
    double ig;
    double igp;
    // scale, usually coincide with pole masses
    double scale;

    bool onlyMasses;
  public:
    MSinput(double mb_, double mW_, double mZ_, double mH_, double mt_) : 
      imb(mb_), imW(mW_), imZ(mZ_), imH(mH_), imt(mt_), onlyMasses(true)
    {
    }
    // With vev
    MSinput(double mb_, double mW_, double mZ_, double mH_, double mt_, double v_, double scale_) : 
      imb(mb_), imW(mW_), imZ(mZ_), imH(mH_), imt(mt_), iv(v_), scale(scale_), onlyMasses(false)
    {
    }
    bool operator < (const MSinput& b) const 
    {
      if(imb < b.mb()) return true;
      else if(imb > b.mb()) return false;
      else
      
        if(imW < b.mW()) return true;
        else if(imW > b.mW()) return false;
        else
        
          if(imZ < b.mZ()) return true;
          else if(imZ > b.mZ()) return false;
          else
        
            if(imH < b.mH()) return true;
            else if(imH > b.mH()) return false;
            else
              if(imt < b.mt()) return true;
              else if(imt > b.mt()) return false;
      return false;
    }

    // Factory to construct from running masses
    static MSinput fromMasses(double mb, double mW, double mZ, double mH, double mt)
    {
      return MSinput(mb, mW, mZ, mH, mt);
    }

    // Factory to construct from running couplings
    static MSinput fromCouplings(double g1,     // U(1)  g1 = sqrt(5/3)*gp
                                 double g2,     // SU(2) g2 = g
                                 double yb, 
                                 double yt, 
                                 double lam, 
                                 double mphi,  //Higgs mass parameter
                                                    //normalized as mphi=Mh at
                                                    //tree level
                                 double scale) // Input scale
                                 
    {
      
      double vev = mphi/sqrt(2.*lam);
      
      double mb = vev*yb/sqrt(2);
      double mW = vev*g2/2.;
      double mZ = sqrt(g2*g2+g1*g1*3./5.)*vev/2.;
      double mH = mphi;
      double mt = vev*yt/sqrt(2);
      return MSinput(mb, mW, mZ, mH, mt, vev, scale);
    }

    // Scale mu^2
    double Q2() const
    {
      if(! onlyMasses)
        return scale;
      else
        throw std::logic_error("ERROR: not availbale when constructed from masses");
    }
    double vev() const
    {
      if(! onlyMasses)
        return iv;
      else
        throw std::logic_error("ERROR: not availbale when constructed from masses");
    }
    // constants 
    double g() const
    {
      return 2.*imW/vev();
    }
    double gp() const
    {
      return 2.*sqrt(imZ*imZ-imW*imW)/vev();
    }
    double alpha() const
    {
      return gp()*gp()*g()*g()/(gp()*gp()+g()*g())/4./Pi;
    }  
    // m^2
    double mmb() const
    {
      return imb*imb;
    }
    double mmW() const
    {
      return imW*imW;
    }
    double mmZ() const
    {
      return imZ*imZ;
    }
    double mmH() const
    {
      return imH*imH;
    }
    double mmt() const
    {
      return imt*imt;
    }

    // Weinberg trigonometric
    double cW() const
    {
      return imW/imZ;
    }
    double ccW() const
    {
      return pow(imW/imZ,2);
    }
    double sW() const
    {
      return sqrt(1-ccW());
    }
    double ssW() const
    {
      return 1-ccW();
    }


    double mb() const
    {
      return imb;
    }
    double mW() const
    {
      return imW;
    }
    double mZ() const
    {
      return imZ;
    }
    double mH() const
    {
      return imH;
    }
    double mt() const
    {
      return imt;
    }

    // modification
    MSinput& setmb(double mb)
    {
      imb = mb;
      return *this;
    }
    MSinput& setmW(double mW)
    {
      imW = mW;
      return *this;
    }
    MSinput& setmZ(double mZ)
    {
      imZ = mZ;
      return *this;
    }
    MSinput& setmH(double mH)
    {
      imH = mH;
      return *this;
    }
    MSinput& setmt(double mt)
    {
      imt = mt;
      return *this;
    }

  };

  typedef MSinput MS;
  typedef OSinput OS;

} // namespace mr

#endif // __SMINPUT_HPP__

