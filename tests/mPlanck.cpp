#include "catch.hpp"
#include "mr.hpp"

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Running up to the Planck scale",
 "[1307.3536][g1][g2][gs][yt][lam]"
 )
{
  
  CouplingsVevMu av(
                    5./3.*pow(0.35830/4./Pi,2), // GUT normalization
                    pow(0.64779/4./Pi,2),
                    pow(1.1666/4./Pi,2),
                    pow(0.93690/4./Pi,2),
                    pow(0.0/4./Pi,2),
                    pow(0.0/4./Pi,2),
                    0.12604*pow(4.*Pi,-2),
                    0*246,
                    // 132.03/sqrt(2.*0.12604),
                    131.55,
                    pow(173.34,2),
                    3);
  
  
  long double mmPlanck = pow(1.2209,2) * pow(10.,2*19); 
  state_type avMpl = av(mmPlanck);
  
  REQUIRE( (4.*Pi*sqrt(avMpl[0])) == Approx( 0.6154 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[1])) == Approx( 0.5055 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[2])) == Approx( 0.4873 ).epsilon( 0.0001 ) );
  REQUIRE( (4.*Pi*sqrt(avMpl[3])) == Approx( 0.3825 ).epsilon( 0.0001 ) );
  REQUIRE( (16.*Pi*Pi*avMpl[6])   == Approx(-0.0143 ).epsilon( 0.0001 ));
  REQUIRE( avMpl[8]               == Approx( 129.4  ).epsilon( 0.0001 ));


  std::cout << "v = " << avMpl[7] << std::endl;
  std::cout << "vc= " << avMpl[8]/sqrt(2.*16.*Pi*Pi*avMpl[6]) << std::endl;

  // state_type av1000 = av(mmPlanck);

}
