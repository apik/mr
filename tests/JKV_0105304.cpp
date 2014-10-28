#include "catch.hpp"
#include "mr.hpp"

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Bosonic part of two-loop EW corrections to the W,Z-boson masses",
 "[0105304][W][Z]"
 )
{
  MSinput mi = MSinput(0., 80.419, 91.188, 125.6, 174.3);
  SECTION( "One-loop" )
    {
      SECTION( "W,\\mu=MW" ) 
        {
          WW<MS> dMW_at_mu_eq_MW  = WW<MS>(mi, mi.mmW());
          
          REQUIRE( dMW_at_mu_eq_MW.m10().real() == Approx( 103.739 ) );
        }
      SECTION( "Z,\\mu=MZ" ) 
        {
          ZZ<MS> dMZ_at_mu_eq_MZ  = ZZ<MS>(mi, mi.mmZ());
          
          REQUIRE( dMZ_at_mu_eq_MZ.m10().real() == Approx( 35.446 ) );
        }

    }
  SECTION( "Two-loop" )
    {
      // Due to expansion in 1/mH we need enormously large mH = 800 GeV
      mi.setmH(800);
      // and small s_W = 0.1
      mi.setmW(mi.mZ()*sqrt(1 - 0.1));

      SECTION( "W,\\mu=MW" ) 
        {
          
          WW<MS> dMW_at_mu_eq_MW  = WW<MS>(mi, mi.mmW());
          
          REQUIRE( dMW_at_mu_eq_MW.m20(0, 0).real() == Approx( 1.95274e6 ) );
        }
      SECTION( "Z,\\mu=MZ" ) 
        {
          ZZ<MS> dMZ_at_mu_eq_MZ  = ZZ<MS>(mi, mi.mmZ());
       
          REQUIRE( dMZ_at_mu_eq_MZ.m20(0, 0).real() == Approx( 1.95905e6 ).epsilon(0.00001) );
        }

    }
}
