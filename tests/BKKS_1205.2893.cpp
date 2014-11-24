#include "catch.hpp"
#include "mr.hpp"

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Mixed QCD EW corrections to the Yukawa-top and Higgs self-coupling",
 "[1205.1893][top][Higgs]"
 )
{
  Approx approx = Approx::custom().epsilon( 0.00001 );
  
  // Important Mb=0, but for comparison 
  // mmd3->0 manualy in deltayta[...]
  OSinput oi = OSinput(0., 80.399, 91.1876, 125.6, 173.5);
    
  SECTION( "top,\\mu=MZ" ) 
    {
      tt dMtZ  = tt(oi, oi.MMZ());
            
      REQUIRE( dMtZ.x10().real() == Approx( -42.1506 ) );
      REQUIRE( dMtZ.x11().real() == Approx(  975.307 ) );

      REQUIRE( dMtZ.y10().real() == Approx( -20.4608 ) );
      REQUIRE( dMtZ.y11().real() == Approx( -167.7   ) );
    }
  SECTION( "top,\\mu=Mt" ) 
    {
      tt dMtt  = tt(oi, oi.MMt());
            
      REQUIRE( dMtt.x10().real() == Approx( 109.666 ) );
      REQUIRE( dMtt.x11().real() == Approx( -417.114 ) );

      REQUIRE( dMtt.y10().real() == Approx( 2.14111 ) );
      REQUIRE( dMtt.y11().real() == Approx( -80.134 ) );
    }
  SECTION( "Higgs,\\mu=Mz" ) 
    {
      HH<OS> dMHZ  = HH<OS>(oi, oi.MMZ());
            
      REQUIRE( dMHZ.x10().real() == Approx( 159.93 ) );
      REQUIRE( dMHZ.x11().real() == Approx( -1527.21 ) );

      REQUIRE( dMHZ.y10().real() == Approx(  203.31 ) );
      REQUIRE( dMHZ.y11().real() == Approx( -3805.1 ) );
    }
  SECTION( "Higgs,\\mu=Mt" ) 
    {
      HH<OS> dMHt  = HH<OS>(oi, oi.MMt());

      REQUIRE( dMHt.x10().real() == Approx( 196.631 ) );
      REQUIRE( dMHt.x11().real() == Approx( -1480.94 ) );

      REQUIRE( dMHt.y10().real() == Approx( -18.4183 ) );
      REQUIRE( dMHt.y11().real() == Approx( -1953.91 ) );
    }
  SECTION( "\\delta-r,\\mu=Mz" ) 
    {
      dr ddrZ  = dr(oi, oi.MMZ());
            
      REQUIRE( ddrZ.dr10().real() == Approx( -43.3796 ) );
      REQUIRE( ddrZ.dr11().real() == Approx(  2277.89 ) );
    }
  SECTION( "\\delta-r,\\mu=Mt" ) 
    {
      dr ddrt  = dr(oi, oi.MMt());
            
      REQUIRE( ddrt.dr10().real() == Approx( 215.049 ) );
      REQUIRE( ddrt.dr11().real() == Approx( 472.97 ) );
    }
}
