#ifndef __SMINPUT_HPP__
#define __SMINPUT_HPP__
class SMinput
{
  long double iMb;
  long double iMW;
  long double iMZ;
  long double iMH;
  long double iMt;

public:
  SMinput(long double Mb_, long double MW_, long double MZ_, long double MH_, long double Mt_) : 
    iMb(Mb_), iMW(MW_), iMZ(MZ_), iMH(MH_), iMt(Mt_)
  {
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
#endif // __SMINPUT_HPP__
