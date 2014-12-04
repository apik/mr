#include <iostream>
#include <Eigen/Dense>
#include "mr.hpp"

using namespace Eigen;

std::complex<long double> det(const long double & a, const long double & b, const long double & c)
{
  return 1./(a*a + b*b + c*c - 2*a*b - 2*b*c - 2*c*a);
}


int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS;

      long double xxx = 1.;

      // Compare with:
      std::vector<OSinput> KV;
      KV.push_back(OSinput(4.4, 80.385, 91.1876, 124, 173.5));
      KV.push_back(OSinput(4.4, 80.385, 91.1876, 125, 173.5));
      KV.push_back(OSinput(4.4, 80.385, 91.1876, 126, 173.5));

      OSinput KVphys = OSinput(4.4, 80.385, 91.1876, 126, 173.5);

      // PDG 2014
      std::vector<OSinput> KPV;
      KPV.push_back(OSinput(4.4, 80.385, 91.1876, 125.7, 173.2));
      KPV.push_back(OSinput(4.4, 80.385, 91.1876, 125.7, 173.2));
      KPV.push_back(OSinput(4.4, 80.385, 91.1876, 125.7, 173.2));

      OSinput KPVphys = OSinput(4.4, 80.385, 91.1876, 125.7, 173.2);


      // tt dMt(KPVphys, 1);
      // yTTmu1gl ttmu1(KPVphys, 1);
      
      
      // long double y = 1.25, x = 2.35, z = 3.45, u = 4.45, v = 5.45, s = 100.555, QQ = 100.;

      
      // std::cout << " B    " << Tsil::B(x, y, s, QQ ) - (Tsil::B(x, y, s, 1. ) + log(QQ)) << std::endl; 
      // std::cout << " Beps " << Tsil::Beps(x,y,s, QQ ) - ( Tsil::Beps(x,y,s, 1. ) + Tsil::B(x,y,s, 1. )*log(QQ) + 1./2.*pow(log(QQ),2)) << std::endl; 


      // std::cout << " Re(B)    " << std::real(Tsil::B(x, y, s, QQ )) - (std::real(Tsil::B(x, y, s, 1. )) + log(QQ)) << std::endl; 
      // std::cout << " Re(Beps) " << std::real(Tsil::Beps(x,y,s, QQ )) - ( std::real(Tsil::Beps(x,y,s, 1. )) + std::real(Tsil::B(x,y,s, 1. ))*log(QQ) + 1./2.*pow(log(QQ),2)) << std::endl; 




      // std::cout << Tsil::A(3., 100. ) - (Tsil::A(3., 1. ) - 3.*log(100.)) << std::endl; 
      // std::cout << Tsil::Aeps(3., 100. ) - ( Tsil::Aeps(3., 1. ) + Tsil::A(3., 1. )*log(100.) - 1./2.*3.*pow(log(100.),2)) << std::endl; 


      // /*

      // id Ifin(x?,y?,z?) = IfinOne(x,y,z)
      //   + (- x - y - z + 2*AtsilOne(x) + 2*AtsilOne(y) + 2*AtsilOne(z) )*Log(mu2/one)
      //   + (- x - y - z)*Log(mu2/one)^2;

      // id Sfin(A(qq?),x?,y?,z?) = SfinOne(A(qq),x,y,z)
      //   + (- x - y - z + qq/2 + 2*AtsilOne(x) + 2*AtsilOne(y) + 2*AtsilOne(z) )*Log(mu2/one)
      //   + (- x - y - z)*Log(mu2/one)^2;


      // id Tfin(A(qq?),x?,y?,z?) = TfinOne(A(qq),x,y,z)
      //   + (-2)*( AtsilOne(x)/x + 1/2 )*Log(mu2/one)
      //   + (+1)*Log(mu2/one)^2;

      // id Ufin(A(qq?),x?,y?,z?,u?) = UfinOne(A(qq),x,y,z,u)
      //   + (+2)*( BtsilOne(A(qq),x,y) + 1/2 )*Log(mu2/one)
      //   + (+1)*Log(mu2/one)^2;

      // id Vfin(A(qq?),x?,y?,z?,u?) = VfinOne(A(qq),x,y,z,u)
      //   + (+2)*( (qq+x-y)*(BtsilOne(A(qq),x,y) - 1) + 2*AtsilOne(x)
      //            + (qq-x-y)*AtsilOne(y)/y )*det(qq,x,y,-1)*Log(mu2/one);

      // id Mfin(?a) = MfinOne(?a); 

      // */


      // std::cout << "I : " << Tsil::I2(1.,2.,3., 100. ) - (
      //                                           Tsil::I2(1.,2.,3., 1. )
      //                                           + (- 1 - 2 - 3 + 2.*Tsil::A(1.,1.) + 2.*Tsil::A(2.,1.) + 2.*Tsil::A(3.,1.) )*log(100.)
      //                                           + (- 1 - 2 - 3)*pow(log(100.),2)
      //                                           ) << std::endl; 
  
      
      // long double mu2 = 100.;
      
      // // Tsil t1(1., 2., 3., 4., 5., 1., 6.);
      // // Tsil tmu(1., 2., 3., 4., 5., mu2, 6.);


      // Tsil t1(x, y, z, u, v, 1., s);
      // Tsil tmu(x, y, z, u, v, mu2, s);
      

      // std::cout << "U : " << tmu.Uxzuv(0) -(
      //                             t1.Uxzuv(0)
      //                             + (+2.)*( Tsil::B(x,z,s,1.) + 1./2. )*log(mu2)
      //                             + (+1.)*pow(log(mu2),2)
      //                             ) << std::endl;



      // /*
      //         id Vfin(A(qq?),x?,y?,z?,u?) = VfinOne(A(qq),x,y,z,u)
      //   + (+2)*( (qq+x-y)*(BtsilOne(A(qq),x,y) - 1) + 2*AtsilOne(x)
      //            + (qq-x-y)*AtsilOne(y)/y )*det(qq,x,y,-1)*Log(mu2/one);

      //  */
      // std::cout << "V : " << tmu.Vxzuv(0) -(
      //                                       t1.Vxzuv(0)
      //                                       + (+2.)*( (s+x-z)*(Tsil::B(x,z,s,1.) - 1.) + 2.*Tsil::A(x,1.)
      //                                                 + (s-x-z)*Tsil::A(z,1.)/z )*det(s,x,z)*log(mu2)
      //                                       ) << std::endl;
      
      // std::cout << "T : " << tmu.Txuv(0) -(
      //                                      t1.Txuv(0)
      //                                      + (-2.)*( Tsil::A(x,1.)/x + 1./2. )*log(mu2)
      //                                      + (+1.)*pow(log(mu2),2)
      //                                      ) << std::endl;

      // std::cout << "S : " << tmu.Suxv(0) -(
      //                                      t1.Suxv(0) 
      //                                      + (- u - x - v + s/2. + 2.*Tsil::A(u,1.) + 2.*Tsil::A(x,1.) + 2.*Tsil::A(v,1.) )*log(mu2)
      //                                      + (- u - x - v)*pow(log(mu2),2)
      //                                      ) << std::endl;
      
      

      return 0;


      // KV.push_back(OSinput(4.89, xxx*80.385, xxx*91.1874, 193.5, 193.5));
      alphaMt  = 0.00781592;
      // alphaMt  = 0.00779305;
      alphaS   = 0.1184;
      
      AlphaS as;
      
      long double alphaMZ = 1./137.035999;

      

      // Btsil(A(mmt), 0, mmW, eps) = Re(Btsil(A(mmt)), 0, mmW, eps)
      //   + 2*IPi*(mmt-mmW)/mmt*( -2 + 2*log(mmt-mmW) - log(mmt) ); 
      
      // for (std::vector<OSinput>::iterator it = KV.begin(); it != KV.end(); ++it)
      //   {
      //     tt dMt  = tt(*it, 2*it->MMt());

      //     dr dDR  = dr(*it, 2*it->MMt());

      //     std::cout << "Test B      : " <<   pow(Tsil::B(0,80.*80.,173.*173.,1.),2) << std::endl;
      //     std::cout << "Test B      : " <<   pow(Tsil::B(0,80.*80.,173.*173.,1.).real(),2) - Pi*Pi*pow((173.*173. - 80.*80.)/173./173.,2) << std::endl;

      //     return 0;
      //     // std::cout << "Test B rhs  : " <<  
      //     //   Re(Tsil::Beps(0,it->MMW(),it->MMt(),it->MMt()))
      //     //   + std::complex<long double>(0,1)*Pi*(it->MMt()-it->MMW())/it->MMt()*( -2 + 2*log((it->MMt()-it->MMW())/it->MMt()) - log(it->MMt()/it->MMt()) ) << std::endl;


      //     std::cout << "Test B      : " <<   Tsil::Beps
      //       (0,it->MMW(),it->MMt(),it->MMt()) << std::endl;
      //     std::cout << "Test B rhs  : " <<  
      //       Re(Tsil::Beps(0,it->MMW(),it->MMt(),it->MMt()))
      //       + std::complex<long double>(0,1)*Pi*(it->MMt()-it->MMW())/it->MMt()*( -2 + 2*log((it->MMt()-it->MMW())/it->MMt()) - log(it->MMt()/it->MMt()) ) << std::endl;
          
          
      //     std::cout << "Mh= " << it->MH()  << std::endl;
      //     std::cout << "as(MMt) = " << as(it->MMt()) << std::endl;          
      //     std::cout << "\t1-loop \\alpha         " << pow(xxx,2)*it->Mt()*alphaMt/4./Pi*dMt.m10() << std::endl;
      //     std::cout << "\t1-loop \\alpha_S       " << it->Mt()*as(it->MMt())/4./Pi*dMt.m01() << std::endl;
      //     std::cout << "\t2-loop \\alpha*\\alpha_S" << pow(xxx,2)*it->Mt()*alphaMt/4./Pi*as(it->MMt())/4./Pi*dMt.m11() << std::endl;
      //     std::cout << "\t2-loop \\alpha^2  g.l. " << pow(xxx,4)*// it->Mt()*
      //       pow(alphaMt/4./Pi,0)*dMt.mgl20() << std::endl;
      //     // std::cout << "\t2-loop \\alpha^2       " << pow(xxx,4)*it->Mt()*pow(alphaMt/4./Pi,2)*dMt.m20() << std::endl;

      //     // std::cout << "\tB2        top     g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),2)*dMt.mgl20() << std::endl;
      //     std::cout << "\tB2        top Y   g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),0)*dMt.mygl20() << std::endl;
      //     // std::cout << "\tB2        top     1 1  " << pow(it->MMW()*it->SSW()/it->MMt(),1)*dMt.m11() << std::endl;
          
      //     std::cout << "Delta-r check:" << std::endl;
      //     std::cout << "\t *****Delta-r   1-loop   top Y   g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),0)*dMt.mygl10() - dMt.mgl10()+1./2.*dDR.drgl10()
      //               << std::endl;
      //     std::cout << "\t *****Delta-r   2-loop   top Y   g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),0)*dMt.mygl20() - dMt.mgl20() + 1./2.*dDR.drgl10()*dMt.mgl10() - 3./8.*pow(dDR.drgl10(),2)+1./2.*dDR.drgl20()
      //               << std::endl;

      //     std::cout << "\n\n\tDelta-r   1-loop   top Y   g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),0)*dDR.drgl10() << std::endl;
      //     std::cout << "\n\n\tDelta-r   2-loop   top Y   g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),0)*dDR.drgl20() << std::endl;
          
      //     // dMt.test();

      //     tt mtA1  = tt(*it, 2*it->MMt());
      //     tt mtA2  = tt(*it, 4*it->MMt());

      //     long double B2 = pow(it->MMW()*it->SSW()/it->MMt(),2)*dMt.mygl20().real();
      //     long double d1 = pow(it->MMW()*it->SSW()/it->MMt(),2)*mtA1.mygl20().real() - B2;
      //     long double d2 = pow(it->MMW()*it->SSW()/it->MMt(),2)*mtA2.mygl20().real() - B2;


      //     Matrix2d A ;
      //     Vector2d b ;
          
      //     A << 
      //       pow(log(1./2.),2), pow(log(1./2.),1),
      //       pow(log(1./4.),2), pow(log(1./4.),1);
          
      //     b << 
      //       d1, d2;
          
      //     FullPivLU< Matrix2d > lu( A ) ;
      //     Vector2d x = lu.solve( b ) ;
          
      //     std::cout << "Solution: " << std::endl <<  "A_2,2 = " << x[0] << std::endl;
      //     std::cout << "          " << std::endl <<  "A_2,1 = " << x[1] << std::endl;
      //     std::cout << "          " << std::endl <<  "B_2   = " << B2   << std::endl;
          
      //     std::cout << "\tB2        top Y   g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),0)*dMt.mgl20() << std::endl;
          

      //     bb dMb  = bb(*it, it->MMt()); // see eq.13
      //     std::cout << "\tB2     bottom          " << pow(it->MMW()*it->SSW()/it->MMt(),1)*dMb.m11() << std::endl;
      //     std::cout << "\tB2     bottom Yukawa   " << pow(it->MMW()*it->SSW()/it->MMt(),1)*dMb.my11() << std::endl;


      //     std::cout << "\tB2     bottom          " << pow(it->MMW()*it->SSW()/it->MMt(),2)*dMb.m20() << std::endl;
      //     std::cout << "\tB2     bottom Yu       " << pow(it->MMW()*it->SSW()/it->MMt(),2)*dMb.my20() << std::endl;


      //     std::cout << "\tB2     bottom     g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),2)*dMb.mgl20() << std::endl;
      //     std::cout << "\tB2     bottom Yu  g.l. " << pow(it->MMW()*it->SSW()/it->MMt(),2)*dMb.mygl20() << std::endl;
          
      //     std::pair<long double,long double>  sumEXECT = dMt.test2(0,0);

      //     std::cout << " EXACT SUM = " << sumEXECT.first*pow(it->MMt(),0)  << " exp = " << sumEXECT.second << std::endl;

      //     std::cout << " EXACT DIFF = " << sumEXECT.first - sumEXECT.second << std::endl;

          
      //     // dMt.test();

      //     std::cout << std::endl << "mas3 test:" << std::endl;
      //     TsilST h00 = TsilST(it->MMH(), 0, 0, it->MMt(), it->MMt());
      //     std::cout << "Txuv  = " << h00.Txuv(0) << std::endl;
      //     std::cout << "Suxv  = " << h00.Suxv(0) << std::endl;

        
      
 
      
      // // I(h,t,t) - CORRECT
      // // I(0,h,t) 
      // std::cout << "I(0,W,t;mmt) = " << Tsil::I2(0,pow(80,2),173*173,pow(173,2)) << std::endl; 
      // std::cout << "I(W,t,t;mmt) = " << Tsil::I2(pow(80,2),pow(173,2),173*173,pow(173,2)) << std::endl; 
      // std::cout << "I(W,Z,t;mmt) = " << Tsil::I2(pow(80,2),pow(91,2),173*173,pow(173,2)) << std::endl; 

      // std::cout << "I(0,W,t;2*mmt) = " << Tsil::I2(0,pow(80,2),173*173,2*pow(173,2)) << std::endl; 
      // std::cout << "I(W,t,t;2*mmt) = " << Tsil::I2(pow(80,2),pow(173,2),173*173,2*pow(173,2)) << std::endl; 
      // std::cout << "I(W,Z,t;2*mmt) = " << Tsil::I2(pow(80,2),pow(91,2),173*173,2*pow(173,2)) << std::endl; 

      // // Masters test

      // Tsil mTest = Tsil(1, 1, 1, 1, 1, 1, 1);

      // std::cout << "M     = " << mTest.M(0) << std::endl;
      // std::cout << "Uzxyv = " << mTest.Uzxyv(0) << std::endl;
      // std::cout << "Uuyxv = " << mTest.Uuyxv(0) << std::endl;
      // std::cout << "Uxzuv = " << mTest.Uxzuv(0) << std::endl;
      // std::cout << "Uyuzv = " << mTest.Uyuzv(0) << std::endl;
      // std::cout << "Tvyz  = " << mTest.Tvyz(0) << std::endl;
      // std::cout << "Tuxv  = " << mTest.Tuxv(0) << std::endl;
      // std::cout << "Tyzv  = " << mTest.Tyzv(0) << std::endl;
      // std::cout << "Txuv  = " << mTest.Txuv(0) << std::endl;
      // std::cout << "Svyz  = " << mTest.Svyz(0) << std::endl;
      // std::cout << "Suxv  = " << mTest.Suxv(0) << std::endl;
      // std::cout << "Vzxyv = " << mTest.Vzxyv(0) << std::endl;
      // std::cout << "Vuyxv = " << mTest.Vuyxv(0) << std::endl;
      // std::cout << "Vxzuv = " << mTest.Vxzuv(0) << std::endl;
      // std::cout << "Vyuzv = " << mTest.Vyuzv(0) << std::endl;
      // std::cout << "I     = " << mTest.I2(1,1,1,1) << std::endl;
      

      // Tsil m00000 = Tsil(0, 0, 0, 0, 0, -1, 1);
      // std::cout << "\t\tM(0,0,0,0,0)     = " << m00000.M(0) << std::endl;

}

catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

