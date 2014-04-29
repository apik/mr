#ifndef __SMINPUT_HPP__
#define __SMINPUT_HPP__
class SMinput
{
  long double MW;
  long double MZ;
  long double MH;
  long double Mt;

public:
  SMinput(long double MW_, long double MZ_, long double MH_, long double Mt_) : 
    MW(MW_), MZ(MZ_), MH(MH_), Mt(Mt_)
  {
  }
  long double MMW() const
  {
    return MW*MW;
  }
  long double MMZ() const
  {
    return MZ*MZ;
  }
  long double MMH() const
  {
    return MH*MH;
  }
  long double MMt() const
  {
    return Mt*Mt;
  }

  long double mw() const
  {
    return MW;
  }
  long double mz() const
  {
    return MZ;
  }
  long double mh() const
  {
    return MH;
  }
  long double mt() const
  {
    return Mt;
  }

};
#endif // __SMINPUT_HPP__
