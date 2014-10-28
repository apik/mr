#include "catch.hpp"
#include "mr.hpp"

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Full two-loop EW corrections to the W,Z-boson masses MW=MZ case",
 "[W][Z]"
 )
{
  MSinput mi = MSinput(0., 80.419, 91.188, 125.6, 174.3);
  // Due to expansion in 1/mH we need enormously large mH = 800 GeV
  mi.setmH(1800);
  // and small s_W = 0.1
  mi.setmW(mi.mZ()*sqrt(1 - 0.001));
  mi.setmt(0.4);
  
  WW<MS> dMW_at_mu_eq_MZ  = WW<MS>(mi, mi.mmZ());
  ZZ<MS> dMZ_at_mu_eq_MZ  = ZZ<MS>(mi, mi.mmZ());
  
  REQUIRE( dMW_at_mu_eq_MZ.m20().real() == Approx( dMZ_at_mu_eq_MZ.m20().real() ).epsilon(0.001) );
}
