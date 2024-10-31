#include "catch.hpp"
#include "mr.hpp"

using namespace mr;

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Strong coupling running and decoupling",
 "[as]"
 )
{
  
  OSinput oi(pdg2014::mb,pdg2014::MW,pdg2014::MZ,pdg2014::MH,pdg2014::Mt);
  AlphaS aS(oi);
  
  double aMt5,aMt6,a100k6;
  
  aMt5 = run(pdg2014::asMZ, pdg2014::MZ, pdg2014::Mt, 5 , 4);
      
  aMt6 = as5nf2as6nf(pdg2014::Mt, pdg2014::Mt, aMt5, /* nl= */5, 3);
  
  a100k6 = run(aMt6, pdg2014::Mt, 100000, 6 );
  
  SECTION( "\\mu=Mt,nf=5" )     // Only running
    {
      OSinput oi(pdg2014::mb,pdg2014::MW,pdg2014::MZ,pdg2014::MH,pdg2014::Mt);
      AlphaS aSnoDec(oi, pdg2014::asMZ, 4, 5);

      REQUIRE( aSnoDec(pdg2014::Mt) == Approx( 0.108002 ) );
    }
  SECTION( "\\mu=Mt,nf=6" )     // Running wth decoupling at Mt
    {
      REQUIRE( aS(pdg2014::Mt) == Approx( 0.108057 ) );
    }
  SECTION( "\\mu=mb,nf=5" ) 
    {
      REQUIRE( aS(pdg2014::mb) == Approx( 0.226525 ) );
    }
  SECTION( "\\mu=100 TeV,nf=6" ) 
    {
      REQUIRE( aS(100000) == Approx( 0.0605896 ) );
    }
}
