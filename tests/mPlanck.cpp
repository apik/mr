#include "catch.hpp"
#include "mr.hpp"

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Running up to the Planck scale",
 "[1307.3536][g1][g2][gs][yt][lam]"
 )
{
  
  CouplingsSM<3,3,3,3,-1,-1,3> av(
                                  5./3.*pow(0.35830/4./Pi,2), // GUT normalization
                                  pow(0.64779/4./Pi,2),
                                  pow(1.1666/4./Pi,2),
                                  pow(0.9369/4./Pi,2),
                                  pow(0.0/4./Pi,2),
                                  pow(0.0/4./Pi,2),
                                  0.12604*pow(4.*Pi,-2),
                                  pow(173.34,2),
                                  3);
  
  
  long double mmPlanck = pow(1.2209,2) * pow(10.,2*19); 
  state_type avMpl = av(mmPlanck);
  
  REQUIRE( (4.*Pi*sqrt(avMpl[0])) == Approx( 0.6154 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[1])) == Approx( 0.5055 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[2])) == Approx( 0.4873 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[3])) == Approx( 0.3825 ).epsilon( 0.0001 ) );
  REQUIRE( (16.*Pi*Pi*avMpl[6])   == Approx(-0.0143 ).epsilon( 0.0001 ));
}
