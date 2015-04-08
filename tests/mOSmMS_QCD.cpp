#include "catch.hpp"
#include "mr.hpp"
#include <Eigen/Dense>

///////////////////////////////////////////////////////////////////////////////
TEST_CASE
(
 "QCD relation between OS and MS masses up to 4-loop",
 "[OS][MS][QCD]"
 )
{
  Approx approx = Approx::custom().epsilon( 0.00001 );
  
  OSinput oi = OSinput(pdg2014::Mb, 80.399, 91.1876, 125.6, 173.5);
    
  SECTION( "MS mass from OS one,\\mu=M" ) 
    {

      // We use b-quark due to nl=2*NL
      bb<OS> dMt(oi, oi.MMb());

      SECTION( "Three-loop zm3 = c0 + c1*nl + c2*nl^2" )
        {
          
          Eigen::Matrix3f mM3;
          Eigen::Vector3f x3NL;
          
          x3NL(0) = dMt.x03(0,1);
          x3NL(1) = dMt.x03(1,1);
          x3NL(2) = dMt.x03(2,1);
          
          // c0 + c1*nl + c2*nl^2
          mM3 <<
            1, 0, 0,
            1, 2, 4,
            1, 4, 16;
          
          Eigen::Vector3f ci3 = pow(4,-3)*mM3.colPivHouseholderQr().solve(x3NL);
          
          // NL^0
          REQUIRE(  ci3(0) == Approx( -198.7068 ) );
          // NL^1
          REQUIRE(  ci3(1) == Approx( 26.9239 ) );
          // NL^2
          REQUIRE(  ci3(2) == Approx( -0.65269 ) );
        }

      SECTION( "Four-loop zm4 = c0 + c1*L + c2*L^2 + c3*L^3 + c4*L^4" )
        {
          
          ////
          /// 4-loop
          // 
          Eigen::Matrix<double, 5, 5> mM4;
          Eigen::Matrix<double, 5, 1> x4L;


          bb<OS> dM1(oi, pow(exp(1),1)*oi.MMb());
          bb<OS> dM2(oi, pow(exp(1),2)*oi.MMb());
          bb<OS> dM3(oi, pow(exp(1),3)*oi.MMb());
          bb<OS> dM4(oi, pow(exp(1),4)*oi.MMb());
          bb<OS> dM5(oi, pow(exp(1),5)*oi.MMb());
          
          x4L(0) = dM1.x04(2,1);
          x4L(1) = dM2.x04(2,1);
          x4L(2) = dM3.x04(2,1);
          x4L(3) = dM4.x04(2,1);
          x4L(4) = dM5.x04(2,1);
          
          
          // c0 + c1*nl + c2*nl^2 + c3*nl^3
          mM4 <<
            1, pow(1,1), pow(1,2), pow(1,3), pow(1,4),
            1, pow(2,1), pow(2,2), pow(2,3), pow(2,4),
            1, pow(3,1), pow(3,2), pow(3,3), pow(3,4),
            1, pow(4,1), pow(4,2), pow(4,3), pow(4,4),
            1, pow(5,1), pow(5,2), pow(5,3), pow(5,4);
          
          Eigen::Matrix<double, 5, 1> ci4 = pow(4,-4)*mM4.colPivHouseholderQr().solve(x4L);
          
          // L^0
          REQUIRE(  ci4(0) == Approx( -1267.0 ) );
          // L^1
          REQUIRE(  ci4(1) == Approx( -500.23 ) );
          // L^2
          REQUIRE(  ci4(2) == Approx( -83.390 ) );
          // L^3
          REQUIRE(  ci4(3) == Approx( -9.9563 ) );
          // L^4
          REQUIRE(  ci4(4) == Approx( -0.514033 ) );
          
        }
    }

  SECTION( "OS mass from MS one,\\mu=M" ) 
    {
      
      
      MSinput mi = MSinput::fromMasses(pdg2014::mb, 1, 2, 3, 4);
      // We use b-quark due to nl=2*NL
      bb<MS> dm(mi, mi.mmb());

      SECTION( "Three-loop cm3 = c0 + c1*nl + c2*nl^2" )
        {
          
          Eigen::Matrix3f Mm3;
          Eigen::Vector3f x3NL;
          
          x3NL(0) = dm.x03(0,1);
          x3NL(1) = dm.x03(1,1);
          x3NL(2) = dm.x03(2,1);
          
          // c0 + c1*nl + c2*nl^2
          Mm3 <<
            1, 0, 0,
            1, 2, 4,
            1, 4, 16;
          
          Eigen::Vector3f ci3 = pow(4,-3)*Mm3.colPivHouseholderQr().solve(x3NL);
          
          // NL^0
          REQUIRE(  ci3(0) == Approx( 190.595 ) );
          // NL^1
          REQUIRE(  ci3(1) == Approx( -26.655 ) );
          // NL^2
          REQUIRE(  ci3(2) == Approx( 0.6527 ).epsilon(0.001) );
        }

      SECTION( "Four-loop cm4 = c0 + c1*L + c2*L^2 + c3*L^3 + c4*L^4" )
        {
          
          ////
          /// 4-loop
          // 
          Eigen::Matrix<double, 5, 5> Mm4;
          Eigen::Matrix<double, 5, 1> x4L;
          

          bb<MS> dm1(mi, pow(exp(1),1)*mi.mmb());
          bb<MS> dm2(mi, pow(exp(1),2)*mi.mmb());
          bb<MS> dm3(mi, pow(exp(1),3)*mi.mmb());
          bb<MS> dm4(mi, pow(exp(1),4)*mi.mmb());
          bb<MS> dm5(mi, pow(exp(1),5)*mi.mmb());
          
          x4L(0) = dm1.x04(2,1);
          x4L(1) = dm2.x04(2,1);
          x4L(2) = dm3.x04(2,1);
          x4L(3) = dm4.x04(2,1);
          x4L(4) = dm5.x04(2,1);
          
          
          // c0 + c1*nl + c2*nl^2 + c3*nl^3
          Mm4 <<
            1, pow(1,1), pow(1,2), pow(1,3), pow(1,4),
            1, pow(2,1), pow(2,2), pow(2,3), pow(2,4),
            1, pow(3,1), pow(3,2), pow(3,3), pow(3,4),
            1, pow(4,1), pow(4,2), pow(4,3), pow(4,4),
            1, pow(5,1), pow(5,2), pow(5,3), pow(5,4);
          
          Eigen::Matrix<double, 5, 1> ci4 = pow(4,-4)*Mm4.colPivHouseholderQr().solve(x4L);
          
          // L^0
          REQUIRE(  ci4(0) == Approx( 1224.0 ).epsilon(0.0001) );
          // L^1
          REQUIRE(  ci4(1) == Approx( 601.98 ).epsilon(0.0001) );
          // L^2
          REQUIRE(  ci4(2) == Approx( 134.10 ).epsilon(0.0001) );
          // L^3
          REQUIRE(  ci4(3) == Approx( 28.846 ).epsilon(0.0001) );
          // L^4
          REQUIRE(  ci4(4) == Approx( 3.9648 ).epsilon(0.0001) );
          
        }
    }
  
}
