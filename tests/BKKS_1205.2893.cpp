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
  
  OSinput oi = OSinput(4.4, 80.399, 91.1876, 125.6, 173.5);
    
  SECTION( "top,\\mu=MZ" ) 
    {
      tt dMtZ  = tt(oi, oi.MMZ());
            
      REQUIRE( dMtZ.m10().real() == Approx( -42.1503 ).epsilon( 0.0001 ) );
      REQUIRE( dMtZ.m11().real() == Approx(  975.307 ).epsilon( 0.0001 ) );

      REQUIRE( dMtZ.my10().real() == Approx( -20.4605 ).epsilon( 0.0001 ) );
      REQUIRE( dMtZ.my11().real() == Approx( -167.7   ).epsilon( 0.0001 ) );
    }
  SECTION( "top,\\mu=Mt" ) 
    {
      tt dMtt  = tt(oi, oi.MMt());
            
      REQUIRE( dMtt.m10().real() == Approx( 109.66 ).epsilon( 0.01 ) );
      REQUIRE( dMtt.m11().real() == Approx( -417.114 ).epsilon( 0.01 ) );

      REQUIRE( dMtt.my10().real() == Approx( 2.13494 ).epsilon( 0.0001 ) );
      REQUIRE( dMtt.my11().real() == Approx( -75.9081   ).epsilon( 0.0001 ) );
    }
  SECTION( "Higgs,\\mu=Mz" ) 
    {
      HH dMHZ  = HH(oi, oi.MMZ());
            
      REQUIRE( dMHZ.my10().real() == Approx(  203.31 ) );
      REQUIRE( dMHZ.my11().real() == Approx( -3805.1 ) );
    }
  SECTION( "Higgs,\\mu=Mt" ) 
    {
      HH dMHt  = HH(oi, oi.MMt());
            
      REQUIRE( dMHt.my10().real() == Approx( -18.4183 ) );
      REQUIRE( dMHt.my11().real() == Approx( -1953.91 ) );
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

  
  INFO( "Here is info" );
  // CAPTURE( dMtZ.m11().real() );
  


}
