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

  typedef OSinputTemplate<long double> OSinput;



  class MSinput
  {
    long double imb;
    long double imW;
    long double imZ;
    long double imH;
    long double imt;
    // only MS
    long double iv;
    long double ig;
    long double igp;
    // scale, usually coincide with pole masses
    long double scale;

    bool onlyMasses;
  public:
    MSinput(long double mb_, long double mW_, long double mZ_, long double mH_, long double mt_) : 
      imb(mb_), imW(mW_), imZ(mZ_), imH(mH_), imt(mt_), onlyMasses(true)
    {
    }
    // With vev
    MSinput(long double mb_, long double mW_, long double mZ_, long double mH_, long double mt_, long double v_, long double scale_) : 
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
    static MSinput fromMasses(long double mb, long double mW, long double mZ, long double mH, long double mt)
    {
      return MSinput(mb, mW, mZ, mH, mt);
    }
  
    // Factory to construct from running couplings
    static MSinput fromCouplings(long double g1,     // U(1)  g1 = sqrt(5/3)*gp
                                 long double g2,     // SU(2) g2 = g
                                 long double yb, 
                                 long double yt, 
                                 long double lam, 
                                 long double mu0,   //Higgs mass parameter
                                                    //normalized as mu0=Mh at
                                                    //tree level
                                 long double scale) // Input scale
                                 
    {
      
      long double vev = mu0/sqrt(2.*lam);
      
      long double mb = vev*yb/sqrt(2);
      long double mW = vev*g2/2.;
      long double mZ = sqrt(g2*g2+g1*g1*3./5.)*vev/2.;
      long double mH = mu0;
      long double mt = vev*yt/sqrt(2);
      return MSinput(mb, mW, mZ, mH, mt, vev, scale);
    }

    // Scale mu^2
    long double Q2() const
    {
      if(! onlyMasses)
        return scale;
      else
        throw std::logic_error("ERROR: not availbale when constructed from masses");
    }
    long double vev() const
    {
      if(! onlyMasses)
        return iv;
      else
        throw std::logic_error("ERROR: not availbale when constructed from masses");
    }
    // constants 
    long double g() const
    {
      return 2.*imW/vev();
    }
    long double gp() const
    {
      return 2.*sqrt(imZ*imZ-imW*imW)/vev();
    }
    long double alpha() const
    {
      return gp()*gp()*g()*g()/(gp()*gp()+g()*g())/4./Pi;
    }  
    // m^2
    long double mmb() const
    {
      return imb*imb;
    }
    long double mmW() const
    {
      return imW*imW;
    }
    long double mmZ() const
    {
      return imZ*imZ;
    }
    long double mmH() const
    {
      return imH*imH;
    }
    long double mmt() const
    {
      return imt*imt;
    }

    // Weinberg trigonometric
    long double cW() const
    {
      return imW/imZ;
    }
    long double ccW() const
    {
      return pow(imW/imZ,2);
    }
    long double sW() const
    {
      return sqrt(1-ccW());
    }
    long double ssW() const
    {
      return 1-ccW();
    }


    long double mb() const
    {
      return imb;
    }
    long double mW() const
    {
      return imW;
    }
    long double mZ() const
    {
      return imZ;
    }
    long double mH() const
    {
      return imH;
    }
    long double mt() const
    {
      return imt;
    }

    // modification
    MSinput& setmb(long double mb)
    {
      imb = mb;
      return *this;
    }
    MSinput& setmW(long double mW)
    {
      imW = mW;
      return *this;
    }
    MSinput& setmZ(long double mZ)
    {
      imZ = mZ;
      return *this;
    }
    MSinput& setmH(long double mH)
    {
      imH = mH;
      return *this;
    }
    MSinput& setmt(long double mt)
    {
      imt = mt;
      return *this;
    }

  };

  typedef MSinput MS;
  typedef OSinput OS;

} // namespace mr

#endif // __SMINPUT_HPP__

