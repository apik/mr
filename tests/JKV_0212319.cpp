#include "catch.hpp"
#include "mr.hpp"

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "Full two-loop EW corrections to the W,Z-boson masses",
 "[0212319][W][Z]"
 )
{
  MSinput mi = MSinput(0., 80.419, 91.188, 800, 174.3);
  OSinput oi = OSinput(0., 80.419, 91.188, 125.6, 0.1);
  SECTION( "Two-loop" )
    {
      // Due to expansion in 1/mH we need enormously large mH = 800 GeV
      mi.setmH(1800);
      // and small s_W = 0.1
      // mi.setmW(mi.mZ()*sqrt(1 - 0.01));
      std::cout << mi.mH() << std::endl;
      SECTION( "W,\\mu=MW" ) 
        {
          
          ww dMW_at_mu_eq_MW  = ww(mi, mi.mmW());

          ww dMZ_at_mu_eq_MZ  = ww(mi, mi.mmZ());
          
          tt dMt  = tt(oi, mi.mmW());
          WW dMW  = WW(oi, mi.mmW());
          ZZ dMZ  = ZZ(oi, mi.mmZ());
          // mi.setmt(0);
          
          std::cout << dMW_at_mu_eq_MW.m10(1,0,0).real() << "   --   " <<  dMW_at_mu_eq_MW.m10(0,1,0).real() << std::endl;
          std::cout << dMW_at_mu_eq_MW.m11(1,0,0).real() << "   --   " <<  dMW_at_mu_eq_MW.m11(0,1,0).real() << std::endl;

          
          std::cout << dMW_at_mu_eq_MW.m20(1,0,1).real() << "   --   " <<  dMZ_at_mu_eq_MZ.m20(1,0,1).real() << std::endl;
          std::cout << dMW_at_mu_eq_MW.m20(0,0,1).real() << "   --   " <<  dMZ_at_mu_eq_MZ.m20(0,0,1).real() << std::endl;


          std::cout << std::setprecision(20);

          std::cout << "\n" << dMt.m10(1,0).real() << "   --   " <<  dMt.m10(0,1).real() << std::endl;
          std::cout << "" << dMt.m11(1,0).real() << "   --   " <<  dMt.m11(0,1).real() << std::endl;
          std::cout << "" << dMt.m20(1,0).real() << "   --   " <<  dMt.m20(0,1).real() << std::endl;

          std::cout << "\n" << dMW.m10(1,0).real() << "   --   " <<  dMW.m10(0,1).real() << std::endl;
          std::cout << "" << dMW.m11(1,0).real() << "   --   " <<  dMW.m11(0,1).real() << std::endl;
          std::cout << "" << dMW.m20(1,0).real() << "   --   " <<  dMW.m20(0,1).real() << std::endl;

          std::cout << "\n" << dMZ.m10(1,0).real() << "   --   " <<  dMZ.m10(0,1).real() << std::endl;
          std::cout << "" << dMZ.m11(1,0).real() << "   --   " <<  dMZ.m11(0,1).real() << std::endl;
          std::cout << "" << dMZ.m20(1,0).real() << "   --   " <<  dMZ.m20(0,1).real() << std::endl;
            
          // dMW.test2();
          
          REQUIRE( dMW_at_mu_eq_MW.m11(1,0,0).real() == Approx( dMW_at_mu_eq_MW.m11(0,1,0).real() ).epsilon(0.00001) );
          REQUIRE( dMW_at_mu_eq_MW.m20().real() == Approx( 1.57207e9 ).epsilon(0.001) );
        }
      SECTION( "Z,\\mu=MZ" ) 
        {
          zz dMZ_at_mu_eq_MZ  = zz(mi, mi.mmZ());
          // mi.setmt(0.1);
          REQUIRE( dMZ_at_mu_eq_MZ.m20(1,0).real() == Approx( 1.95306e+09 ).epsilon(0.00001) );

          REQUIRE( dMZ_at_mu_eq_MZ.m20(1,0).real() == Approx( dMZ_at_mu_eq_MZ.m20(0,1).real() ).epsilon(0.00001) );
          REQUIRE( dMZ_at_mu_eq_MZ.m20(1,0).real() == Approx( 1.95306e+09 ).epsilon(0.00001) );
        }

    }
}
