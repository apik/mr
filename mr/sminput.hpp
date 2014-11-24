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
class OSinput
{
  long double iMb;
  long double iMW;
  long double iMZ;
  long double iMH;
  long double iMt;

public:
  OSinput(long double Mb_, long double MW_, long double MZ_, long double MH_, long double Mt_) : 
    iMb(Mb_), iMW(MW_), iMZ(MZ_), iMH(MH_), iMt(Mt_)
  {
  }

  OSinput() 
  {
  }
  
  bool operator < (const OSinput& b) const 
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

  long double MMb() const
  {
    return iMb*iMb;
  }
  long double MMW() const
  {
    return iMW*iMW;
  }
  long double MMZ() const
  {
    return iMZ*iMZ;
  }
  long double MMH() const
  {
    return iMH*iMH;
  }
  long double MMt() const
  {
    return iMt*iMt;
  }

  // Weinberg trigonometric
  long double CW() const
  {
    return iMW/iMZ;
  }
  long double CCW() const
  {
    return pow(iMW/iMZ,2);
  }
  long double SW() const
  {
    return sqrt(1-CCW());
  }
  long double SSW() const
  {
    return 1-CCW();
  }


  long double Mb() const
  {
    return iMb;
  }
  long double MW() const
  {
    return iMW;
  }
  long double MZ() const
  {
    return iMZ;
  }
  long double MH() const
  {
    return iMH;
  }
  long double Mt() const
  {
    return iMt;
  }

  // Modification
  void setMb(long double mb)
  {
    iMb = mb;
  }
  void setMW(long double mW)
  {
    iMW = mW;
  }
  void setMZ(long double mZ)
  {
    iMZ = mZ;
  }
  void setMH(long double mH)
  {
    iMH = mH;
  }
  void setMt(long double mt)
  {
    iMt = mt;
  }

};

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

public:
  MSinput(long double mb_, long double mW_, long double mZ_, long double mH_, long double mt_) : 
    imb(mb_), imW(mW_), imZ(mZ_), imH(mH_), imt(mt_)
  {
  }
  // With vev
  MSinput(long double mb_, long double mW_, long double mZ_, long double mH_, long double mt_, long double v_) : 
    imb(mb_), imW(mW_), imZ(mZ_), imH(mH_), imt(mt_), iv(v_)
  {
  }

  // Factory
  static MSinput fromMasses(long double mb, long double mW, long double mZ, long double mH, long double mt)
  {
    return MSinput(mb, mW, mZ, mH, mt);
  }
  static MSinput fromConsts(long double lam, 
                            long double v, 
                            long double yb, 
                            long double yt, 
                            long double g, // SU(2) 
                            long double gp // U(1)
                            )
  {
    long double mb = v*yb/sqrt(2);
    long double mW = v*g/2.;
    long double mZ = sqrt(g*g+gp*gp)*v/2.;
    long double mH = v*sqrt(2.*lam);
    long double mt = v*yt/sqrt(2);
    return MSinput(mb, mW, mZ, mH, mt, v);
  }
  // constants 
  long double g() const
  {
    return 2.*imW/iv;
  }
  long double gp() const
  {
    return 2.*sqrt(imZ*imZ-imW*imW)/iv;
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
  void setmb(long double mb)
  {
    imb = mb;
  }
  void setmW(long double mW)
  {
    imW = mW;
  }
  void setmZ(long double mZ)
  {
    imZ = mZ;
  }
  MSinput& setmH(long double mH)
  {
    imH = mH;
    return *this;
  }
  void setmt(long double mt)
  {
    imt = mt;
  }

};

typedef MSinput MS;
typedef OSinput OS;
#endif // __SMINPUT_HPP__

