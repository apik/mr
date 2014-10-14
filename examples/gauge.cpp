#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaSMt;

      // Compare with:
      SMinput KVPhys(4.40, 80.385, 91.1876, 125.6, 173.5);
      

      long double alphaSMZ = 0.1184;
      // \mu = Mt
      alphaMt  = 0.00779305;
      // From Degrassi
      alphaMt  = 0.00780216;
      alphaSMt   = 0.1079;

      
      // \mu = Mb
      long double alphaMb  = 0.00784257;
      long double alphaSMb   = 0.1905;

      

      AlphaS as;
      
      long double alphaMZ = 1./127.944;


      // 
      // 
      // Test relation betwee alpha and G Fermi
      // 
      // 
      
      
      WW dW  = WW(KVPhys, KVPhys.MMt());
      ZZ dZ  = ZZ(KVPhys, KVPhys.MMt());

      long double alpha    = alphaMZ;
      long double alphaS   = alphaSMZ;

      long double aHat,Gf;
      
      Gf = 1.16637e-5;

      std::complex<long double> dMyW,dMyZ;
      dMyW = 1;
      dMyZ = 1;

      std::cout << " Input:  1/" << 1./alpha << std::endl;

      // Tree level
      aHat = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);
      std::cout << " Tree:   1/" << 1./aHat << std::endl;

      // 1-loop level
      dMyW += alpha/4./Pi*dW.my10();
      dMyZ += alpha/4./Pi*dZ.my10();

      aHat = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);
      std::cout << " 1-loop: 1/" << 1./aHat << std::endl;

      // 2-loop level EW QCD
      dMyW += alpha/4./Pi*alphaS/4./Pi*dW.my11();
      dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ.my11();

      aHat = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);
      std::cout << " EW*QCD: 1/" << 1./aHat << std::endl;

      // 2-loop level EW
      dMyW += pow(alpha/4./Pi,2)*dW.my20();
      dMyZ += pow(alpha/4./Pi,2)*dZ.my20();

      aHat = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);
      std::cout << " EW*EW:  1/" << 1./aHat << std::endl;

      // // Yukawa W
      // WW dMW  = WW(KVPhys, KVPhys.MMt());
      // std::cout << "[ W ]" << std::endl;
      // std::cout << "Mh= " << KVPhys.MH()  << std::endl;
      // std::cout << "as(MMt) = " << as(KVPhys.MMt()) << std::endl;          
      // std::cout << "\t1-loop \\alpha         " << alphaMZ/4./Pi*dMW.my10() << std::endl;
      // // std::cout << "\t1-loop \\alpha_S       " << alphaSMZ/4./Pi*dMW.my01() << std::endl;
      // std::cout << "\t2-loop \\alpha*\\alpha_S" << alphaMZ/4./Pi*alphaSMZ/4./Pi*dMW.my11() << std::endl;
      // std::cout << "\t2-loop \\alpha^2       " << pow(alphaMZ/4./Pi,2)*dMW.my20() << std::endl;

      // long double toLam =1.+0*KVPhys.MMH()*0.0000116637/sqrt(2);


      // SMinput degr(4.40, 80.384, 91.1876, 125.66, 173.10);

      // HH dMH  = HH(degr, degr.MMt());
      // std::cout << "[ H ]" << std::endl;
      // std::cout << "Mh= " << degr.MH()  << std::endl;
      // std::cout << "as(MMt) = " << as(degr.MMt()) << std::endl;          
      // std::cout << "\t1-loop \\alpha         " << toLam*alphaMt/4./Pi*dMH.lam10() << std::endl;
      // // std::cout << "\t1-loop \\alpha_S       " << alphaSMZ/4./Pi*dMW.my01() << std::endl;
      // std::cout << "\t2-loop \\alpha*\\alpha_S" << toLam*alphaMt/4./Pi*alphaSMt/4./Pi*dMH.lam11() << std::endl;
      // std::cout << "\t2-loop \\alpha^2       " << toLam*pow(alphaMt/4./Pi,2)*dMH.lam20() << std::endl;


      // Yukawa bottom
   //    bb dMb  = bb(degr, degr.MMb());
   //    std::cout << "[ Bottom quark ]" << std::endl;
   //    std::cout << "\t1-loop \\alpha         " << alphaMb/4./Pi*dMb.my10() << std::endl;
   //    std::cout << "\t1-loop \\alpha_S       " << alphaSMb/4./Pi*dMb.my01() << std::endl;
   //    std::cout << "\t2-loop \\alpha*\\alpha_S" << alphaMb/4./Pi*alphaSMt/4./Pi*dMb.my11() << std::endl;
   //    std::cout << "\t2-loop \\alpha^2       " << pow(alphaMb/4./Pi,2)*dMb.my20() << std::endl;
      
      
   //    std::cout << "[ ratios ]" << std::endl;
   //    std::cout << "<10> yt/yb = " << dMb.my10()/dMt.my10() << " mt/mb = " << dMb.m10()/dMt.m10() << std::endl;
   //    std::cout << "<01> yt/yb = " << dMb.my01()/dMt.my01() << " mt/mb = " << dMb.m01()/dMt.m01() << std::endl;
   //    std::cout << "<11> yt/yb = " << dMb.my11()/dMt.my11() << " mt/mb = " << dMb.m11()/dMt.m11() << std::endl;
   //    std::cout << "<20> yt/yb = " << dMb.my20()/dMt.my20() << " mt/mb = " << dMb.m20()/dMt.m20() << std::endl;
      
      

   //    // Test Jegerlehner input
   //    // using 1-loop matching

   //    SMinput inFJ(0,80.385,91.1876,125.5,173.5);

   //    tt topFJ(inFJ, inFJ.MMt());

   //    std::cout << "\n\n \t Jegerlehner input:" << std::endl;
      
   //    std::complex<long double> ytFJ = inFJ.Mt()*(1 
   //                                                + alphaMt/4./Pi*topFJ.m10() 
   //                                                + alphaSMt/4./Pi*topFJ.m01()
   //                                                + alphaMt/4./Pi*alphaSMt/4./Pi*topFJ.m11()
   //                                                )*sqrt(2*sqrt(2)*1.16637e-5);

   //    std::cout << "yT(mt) = " << ytFJ << std::endl;


   //    // Plot Yukawa top
   // Plot1 plotYt("yt", "Yukawa Top", "mH", "\\sigma_\\alpha*\\alpha_S", "a*a_S");
   // long double mHstep  = 20; // GeV
   // long double mHstart = 80; // GeV
   
   // for (int mHi = 0; mHi < 13; mHi++)
   //   {
   //     SMinput DS2l(4.40, 80.385, 91.1876, mHstart + mHi*mHstep, 173.5);
   //     tt dtY  = tt(DS2l, DS2l.MMt());          
       
   //     plotYt.add(DS2l.MH(),alphaMt/4./Pi*alphaSMt/4./Pi*dtY.my11().real());
   //   }


   // // Plot Yukawa bottom
   // Plot1 plotYb("yb", "Yukawa Bottom", "mH", "\\sigma_\\alpha*\\alpha_S", "a*a_S");
   // mHstep  = 10; // GeV
   // mHstart = 1; // GeV
   
   // for (int mHi = 0; mHi < 100; mHi++)
   //   {
   //     SMinput DS2l(4.40, 80.385, 91.1876, 125.6, 173.5);
   //     bb dbY  = bb(DS2l, mHstart + mHi*mHstep);          
       
   //     plotYb.add(mHstart + mHi*mHstep,alphaMb/4./Pi*alphaSMb/4./Pi*dbY.my11().real());
   //   }

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

