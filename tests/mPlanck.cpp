#include "catch.hpp"
#include "mr.hpp"

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Running up to the Planck scale v4",
 "[1307.3536v4][g1][g2][gs][yt][lam]"
 )
{
  
  Couplings<3,3,3,
            3,3,-1,             // No dependence on Yukawa-tau 
            3,3,0>              // And we do not need running VEV
    av(
       5./3.*pow(0.35830/4./Pi,2), // GUT normalization
       pow(0.64779/4./Pi,2),
       pow(1.1666/4./Pi,2),
       pow(0.93690/4./Pi,2),
       pow(0.0/4./Pi,2),
       pow(0.0/4./Pi,2),
       0.12604*pow(4.*Pi,-2),
       131.55,
       0*246,
       pow(173.34,2),
       3);
  
  
  long double mmPlanck = pow(1.2209,2) * pow(10.,2*19); 
  state_type avMpl = av(mmPlanck);
  
  REQUIRE( (4.*Pi*sqrt(avMpl[0])) == Approx( 0.6154 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[1])) == Approx( 0.5055 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[2])) == Approx( 0.4873 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[3])) == Approx( 0.3825 ).epsilon( 0.0001 ) );
  REQUIRE( (16.*Pi*Pi*avMpl[6])   == Approx(-0.0143 ).epsilon( 0.0001 ));
  REQUIRE( avMpl[7]               == Approx( 139.4  ).epsilon( 0.0001 )); // In
                                                                          // original
                                                                          // paper
                                                                          // there
                                                                          // is
                                                                          // mistake 129.4

}


TEST_CASE
(
 "Running up to the Planck scale v2",
 "[1307.3536v2][g1][g2][gs][yt][lam]"
 )
{
  
  Couplings<3,3,3,
            3,3,-1,             // No dependence on Yukawa-tau 
            3,3,0>              // And we do not need running VEV
    av(
       5./3.*pow(0.35761/4./Pi,2), // GUT normalization
       pow(0.64822/4./Pi,2),
       pow(1.1666/4./Pi,2),
       pow(0.93558/4./Pi,2),
       pow(0.0/4./Pi,2),
       pow(0.0/4./Pi,2),
       0.12711*pow(4.*Pi,-2),
       132.03,
       0*246,
       pow(173.10,2),
       3);
  
  
  long double mmPlanck = pow(1.2209,2) * pow(10.,2*19); 
  state_type avMpl = av(mmPlanck);
  
  REQUIRE( (4.*Pi*sqrt(avMpl[0])) == Approx( 0.6133 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[1])) == Approx( 0.5057 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[2])) == Approx( 0.4873 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[3])) == Approx( 0.3813 ).epsilon( 0.0001 ) );
  REQUIRE( (16.*Pi*Pi*avMpl[6])   == Approx(-0.0113 ).epsilon( 0.0001 ));
  REQUIRE( avMpl[7]               == Approx( 140.3  ).epsilon( 0.0001 ));

}
