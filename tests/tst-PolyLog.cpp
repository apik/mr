#include "catch.hpp"
#include "mr.hpp"

using namespace mr;

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Check DiLog and TriLog evaluation implemented in TSIL",
 "[PolyLog]"
 )
{
  
  
  SECTION( "Real arguments" )
    {

      std::complex<double> li2 = mr::Li2( std::complex<double>(1,0) );
      REQUIRE( li2.real() == Approx( 1.64493 ) );
      REQUIRE( li2.imag() == Approx( 0 ) );

      std::complex<double> li3 = mr::Li3( std::complex<double>(1,0) );
      REQUIRE( li3.real() == Approx( 1.20206 ) );
      REQUIRE( li3.imag() == Approx( 0 ) );

    }

  SECTION( "Imaginary arguments" )
    {

      std::complex<double> li2 = mr::Li2( std::complex<double>(0,1) );
      REQUIRE( li2.real() == Approx( -0.205617 ) );
      REQUIRE( li2.imag() == Approx( 0.915966 ) );

      std::complex<double> li3 = mr::Li3( std::complex<double>(0,1) );
      REQUIRE( li3.real() == Approx( -0.112693 ) );
      REQUIRE( li3.imag() == Approx( 0.968946 ) );

    }

  SECTION( "Complex arguments" )
    {

      std::complex<double> li2 = mr::Li2( std::complex<double>(1,1) );
      REQUIRE( li2.real() == Approx( 0.61685 ) );
      REQUIRE( li2.imag() == Approx( 1.46036 ) );

      std::complex<double> li3 = mr::Li3( std::complex<double>(1,1) );
      REQUIRE( li3.real() == Approx( 0.871159 ) );
      REQUIRE( li3.imag() == Approx( 1.26708 ) );

    }


}
