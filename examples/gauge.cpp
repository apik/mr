#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);


      long double x = 0.1;
      TsilST Tt00 = TsilST(x,1.,2.,3.,4.);
      std::cout << "x = " << x << "   T(x,0,0,1) = " << Tt00.Txuv(0) << std::endl;

      return 0;


      long double MMt,MMW,MMZ,MMH,alphaMt,alphaSMt;

      OSinput KVPhys(4.40, 80.385, 91.1875, 125.6, 172.3);
      
      // long double alpha = 1./127.944;
      long double alpha = 1./128.342;
      // long double alpha = 1./127.44;

      long double alphaS = 0.1184;

      WW<OS> dW    = WW<OS>(KVPhys, KVPhys.MMZ());
      ZZ<OS> dZ    = ZZ<OS>(KVPhys, KVPhys.MMZ());
      dr drOS  = dr(KVPhys, KVPhys.MMZ());

      long double r10,r11,r20,w10,w11,w20,z10,z11,z20;

      
      r10 = real(alpha/4./Pi*drOS.dr10());
      r11 = real(alpha/4./Pi*alphaS/4./Pi*drOS.dr11());
      r20 = real(pow(alpha/4./Pi,2)*drOS.dr20());
      
      w10 = real(alpha/4./Pi*dW.m10());
      w11 = real(alpha/4./Pi*alphaS/4./Pi*dW.m11());
      w20 = real(pow(alpha/4./Pi,2)*dW.m20());
      
      z10 = real(alpha/4./Pi*dZ.m10());
      z11 = real(alpha/4./Pi*alphaS/4./Pi*dZ.m11());
      z20 = real(pow(alpha/4./Pi,2)*dZ.m20());
      
      long double rCZ10,rCZ11,rCZ20;

      long double s = KVPhys.SW();
      rCZ10 = r10 - w10 + (pow(s,-2) - 1)*(w10 - z10);
      
      rCZ11 = r11 - w11 + (pow(s,-2) - 1)*(w11 - z11);

      rCZ20 = r20 - w20 + (1/pow(s,2) - 1)*(w20 - z20)

        +(3 - 3*pow(s,-2) + pow(s,-4))*w10*w10 + (pow(s,-4) - pow(s,-2))*z10*z10 + (-2 + pow(s,-2))*r10*w10

        +(-2 + 4*pow(s,-2) - 2*pow(s,-4))*w10*z10 + (1 - pow(s,-2))*r10*z10;

      std::cout << "dr = 1" << std::endl
                << "   + (" << rCZ10 << ")" <<  std::endl
                << "   + (" << rCZ11 << ")" <<  std::endl
                << "   + (" << rCZ20 << ")" <<  std::endl
                << "   = " << 10000*(rCZ10 + rCZ11 + rCZ20 + 7e-4) <<  std::endl;

      std::cout << "GF = " << alpha*Pi/sqrt(2)*(1 + rCZ10 + rCZ11 + rCZ20 + 7e-4)/KVPhys.MMW()/KVPhys.SSW()<< std::endl;
                // << "   + (" << rCZ10 << ")" <<  std::endl
                // << "   + (" << rCZ11 << ")" <<  std::endl
                // << "   + (" << rCZ20 << ")" <<  std::endl
                // << "   = " << 1 + rCZ10 + rCZ11 + rCZ20 <<  std::endl;


      
      return 0;
      // 2-loop level EW

      std::complex<long double> dMW,dMZ,dR;
      dMW = 0;
      dMZ = 0;
      dR  = 0;

      dMW += alpha/4./Pi*dW.m10();
      dMZ += alpha/4./Pi*dZ.m10();
      dR   += alpha/4./Pi*drOS.dr10();
      
      
      // dMW += alpha/4./Pi*alphaS/4./Pi*dW.m11();
      // dMZ += alpha/4./Pi*alphaS/4./Pi*dZ.m11();
      // dR   += alpha/4./Pi*alphaS/4./Pi*drOS.dr11();
      
      // dMW += pow(alpha/4./Pi,2)*dW.m20();
      // dMZ += pow(alpha/4./Pi,2)*dZ.m20();
      // dR  += pow(alpha/4./Pi,2)*drOS.dr20();

      
      long double drCzakon = (dR*(1 - KVPhys.MMW()/KVPhys.MMZ())/dMW/(1 - KVPhys.MMW()/KVPhys.MMZ()*dMW/dMZ) - 1.).real();

      long double drCzakonExp = real(dR + (dMW )*(1./KVPhys.SSW() - 2.) + (dMZ )*(1-1./KVPhys.SSW()));  


      std::cout << "rCZ = " << (drCzakonExp+ 0.05907)*10000. << std::endl;
      std::cout << "SW  = " <<       KVPhys.SSW() << std::endl;

      std::cout << "W  = " <<       dMW << std::endl;
      std::cout << "Z  = " <<       dMZ << std::endl;
      std::cout << "R  = " <<       dR << std::endl;


      // long double x = 0.1;
      // TsilST Tt00 = TsilST(x,0.,0.,1.,1.);
      // std::cout << "x = " << x << "   T(x,0,0,1) = " << Tt00.Txuv(0) << std::endl;
      return 0;
      // Compare with:
      // OSinput KVPhys(4.40, 80.385, 91.1876, 125.6, 173.5);
      

      // long double alphaTree = 1./137.234;

      // long double alphaSMZ = 0.1184;
      // // \mu = Mt
      // alphaMt  = 0.00779305;
      // // From Degrassi
      // alphaMt  = 0.00780216;
      // alphaSMt   = 0.1079;

      
      // // \mu = Mb
      // long double alphaMb  = 0.00784257;
      // long double alphaSMb   = 0.1905;

      


      // alphaGF aGF  = alphaGF(KVPhys, KVPhys.MMZ());
      
      // long double Gf = 1.16637e-5;

      // std::complex<long double> dMyW,dMyZ;
      // dMyW = 1;
      // dMyZ = 1;


      // // Tree level
      // long double alF = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);

      // WW dW    = WW(KVPhys, KVPhys.MMZ());
      // ZZ dZ    = ZZ(KVPhys, KVPhys.MMZ());


      // long double ali[3];
      // ali[0] = alF;
      // ali[1] = alF*(1 + alF/4./Pi*aGF.a10().real());
      // ali[2] = alF*(1 + alF/4./Pi*aGF.a10().real() + alF/4./Pi*alphaSMZ/4./Pi*aGF.a11().real());
      // ali[3] = alF*(1 + alF/4./Pi*aGF.a10().real() + alF/4./Pi*alphaSMZ/4./Pi*aGF.a11().real() + pow(alF/4./Pi,2)*aGF.a20().real());

      // std::cout << std::setprecision(8);
      // std::cout << "1/alpha = " << 1./ali[0]  << " : born" << std::endl      
      //           << "          " << 1./ali[1] - 1./ali[0] << " : EW " << std::endl
      //           << "          " << 1./ali[2] - 1./ali[1] << " : EW * QCD " << std::endl
      //           << "          " << 1./ali[3] - 1./ali[2] << " : EW * EW " << std::endl
      //           << "          " << std::endl
      //           << "          " << 1./ali[3] << " : Total " << std::endl;
      
      
      // std::cout << aGF.a10() << std::endl;
      // std::cout << dW.my10() << std::endl;
      // std::cout << dZ.my10() << std::endl;
      // return 0;

      
      // // long double alphaMZ = 1./127.944;
      // // long double alphaMZ = 1./127.773;
      // long double alphaMZ = 1./126.654;

      // // mZ with mH->0
      // // std::ofstream of("mh0.dat");
      // // for(long double log10mH = 3; log10mH > -7.; log10mH -= 1.)
      // //   {

          
      // //     long double mH = powf(10.,log10mH);
      // //     OSinput inH0(4.40, 80.385, 91.1876,mH , 173.5);
      // //     ZZ dZ    = ZZ(inH0, KVPhys.MMt());
      // //     WW dW    = WW(inH0, KVPhys.MMt());
      // //     // of        << mH << ", " << dZ.m20(0,0).real() << std::endl;
      // //     std::cout << "Log_10(mH) = " << log10mH <<" Mh = " << mH << ", MW = " << dW.m20(0,0) << ", MZ = " << dZ.m20(0,0) << std::endl;
      //     MMH = inH0.MMH();
      //     MMZ = inH0.MMZ();
      //     MMW = inH0.MMW();
      //     long double mu2 = inH0.MMZ();
      //     std::cout << " Bubbles: \n" 
      //               << Tsil::I2(MMH,MMZ,MMZ,mu2) << std::endl
      //               <<Tsil::I2(MMH,MMW,MMW,mu2)  << std::endl
      //               << Tsil::I2(0,MMZ,MMH,mu2)   << std::endl
      //               << Tsil::I2(0,MMW,MMH,mu2)   << std::endl
      //               << Tsil::I2(0,0,MMZ,mu2)     << std::endl
      //               << Tsil::I2(0,0,MMW,mu2)     << std::endl;

      //   }
      // return 0;
      // 
      // 
      // Test relation betwee alpha and G Fermi
      // 
      // 
      
      
      // WW dW    = WW(KVPhys, KVPhys.MMt());
      // ZZ dZ    = ZZ(KVPhys, KVPhys.MMt());
      // tt dt    = tt(KVPhys, KVPhys.MMt());
      // dr drOS  = dr(KVPhys, KVPhys.MMt());




      // WW dW    = WW(KVPhys, pow(1000,2));
      // ZZ dZ    = ZZ(KVPhys, pow(1000,2));
      // tt dt    = tt(KVPhys, pow(1000,2));
      // dr drOS  = dr(KVPhys, pow(1000,2));

      // long double alpha    = alphaTree;
      // long double alphaS   = alphaSMZ;

      // long double aHat,Gf,ait[3];
      
      // Gf = 1.16637e-5;

      // std::complex<long double> dMyW,dMyZ,dR;
      // dMyW = 1;
      // dMyZ = 1;

      // // std::cout << " Input:  1/" << 1./alpha << std::endl;

      // // Tree level
      // alpha = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);

      // ait[0] = aHat;            // First iteration
      // dR   = 1;      
      // std::cout << "\n Tree:   \\alpha   = 1/" << 1./alpha << std::endl;
      // std::cout << " Tree:   \\delta-r = " << dR << std::endl;

      // // Initial input !!!
      // // alpha = 1./137.03599;
      // alpha = alphaMZ;
      // std::cout << "\n \\alpha   = 1/" << 1./alpha << std::endl;


      // // 1-loop level
      // dMyW += alpha/4./Pi*dW.my10();
      // dMyZ += alpha/4./Pi*dZ.my10();
      // dR   += alpha/4./Pi*drOS.dr10();

      // alpha = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);

      // std::cout << "\n 1-loop: 1/" << 1./alpha << std::endl;
      // std::cout << " \\delta-r = " << dR << std::endl;

      // std::cout << " dMW = " << dW.m10().real() << std::endl;
      // std::cout << " dMZ = " << dZ.m10().real() << std::endl;
      // std::cout << " dMt = " << dt.m10().real() << std::endl;

      
      // std::cout << " dMyW = " << dW.my10().real() << std::endl;
      // std::cout << " dMyZ = " << dZ.my10().real() << std::endl;

      

      // // 2-loop level EW QCD
      // dMyW = 1;
      // dMyZ = 1;

      // dMyW += alpha/4./Pi*dW.my10();
      // dMyZ += alpha/4./Pi*dZ.my10();
      // dR   += alpha/4./Pi*drOS.dr10();

      // dMyW += alpha/4./Pi*alphaS/4./Pi*dW.my11();
      // dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ.my11();
      // dR   += alpha/4./Pi*alphaS/4./Pi*drOS.dr11();
      
      
      // alpha = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);
      // std::cout << "\n EW*QCD: 1/" << 1./alpha << std::endl;
      // std::cout << " \\delta-r = " << dR << std::endl;

      // std::cout << " dMW = " << dW.m11().real() << std::endl;
      // std::cout << " dMZ = " << dZ.m11().real() << std::endl;
      // std::cout << " dMyW = " <<dW.my11().real() << std::endl;
      // std::cout << " dMyZ = " <<dZ.my11().real() << std::endl;

      // std::cout << " dMt = " << dt.m11().real() << std::endl;
      // std::cout << " dMyt = " <<dt.my11().real() << std::endl;
      // // 2-loop level EW
      // dMyW = 1;
      // dMyZ = 1;
      
      // dMyW += alpha/4./Pi*dW.my10();
      // dMyZ += alpha/4./Pi*dZ.my10();
      // dR   += alpha/4./Pi*drOS.dr10();
      
      // dMyW += alpha/4./Pi*alphaS/4./Pi*dW.my11();
      // dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ.my11();
      // dR   += alpha/4./Pi*alphaS/4./Pi*drOS.dr11();
      
      // dMyW += pow(alpha/4./Pi,2)*dW.my20();
      // dMyZ += pow(alpha/4./Pi,2)*dZ.my20();
      // dR   += pow(alpha/4./Pi,2)*drOS.dr20();

      // alpha = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);
      // std::cout << "\n EW*EW:  1/" << 1./alpha << std::endl;
      // std::cout << " \\delta-r = " << dR << std::endl;

      // std::cout << " dMW = " << dW.m20().real() << std::endl;
      // std::cout << " dMZ = " << dZ.m20().real() << std::endl;

      // std::cout << " dMyW = " << dW.my20().real() << std::endl;
      // std::cout << " dMyZ = " << dZ.my20().real() << std::endl;

      // std::cout << " dMt = " << dt.m20().real() << std::endl;
      // std::cout << " dMyt = " << dt.my20().real() << std::endl;


      

      


      // // // Yukawa W
      // // WW dMW  = WW(KVPhys, KVPhys.MMt());
      // // std::cout << "[ W ]" << std::endl;
      // // std::cout << "Mh= " << KVPhys.MH()  << std::endl;
      // // std::cout << "as(MMt) = " << as(KVPhys.MMt()) << std::endl;          
      // // std::cout << "\t1-loop \\alpha         " << alphaMZ/4./Pi*dMW.my10() << std::endl;
      // // // std::cout << "\t1-loop \\alpha_S       " << alphaSMZ/4./Pi*dMW.my01() << std::endl;
      // // std::cout << "\t2-loop \\alpha*\\alpha_S" << alphaMZ/4./Pi*alphaSMZ/4./Pi*dMW.my11() << std::endl;
      // // std::cout << "\t2-loop \\alpha^2       " << pow(alphaMZ/4./Pi,2)*dMW.my20() << std::endl;

      // // long double toLam =1.+0*KVPhys.MMH()*0.0000116637/sqrt(2);


      // // OSinput degr(4.40, 80.384, 91.1876, 125.66, 173.10);

      // // HH dMH  = HH(degr, degr.MMt());
      // // std::cout << "[ H ]" << std::endl;
      // // std::cout << "Mh= " << degr.MH()  << std::endl;
      // // std::cout << "as(MMt) = " << as(degr.MMt()) << std::endl;          
      // // std::cout << "\t1-loop \\alpha         " << toLam*alphaMt/4./Pi*dMH.lam
      // 10() << std::endl;
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

   //    OSinput inFJ(0,80.385,91.1876,125.5,173.5);

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
   //     OSinput DS2l(4.40, 80.385, 91.1876, mHstart + mHi*mHstep, 173.5);
   //     tt dtY  = tt(DS2l, DS2l.MMt());          
       
   //     plotYt.add(DS2l.MH(),alphaMt/4./Pi*alphaSMt/4./Pi*dtY.my11().real());
   //   }


   // // Plot Yukawa bottom
   // Plot1 plotYb("yb", "Yukawa Bottom", "mH", "\\sigma_\\alpha*\\alpha_S", "a*a_S");
   // mHstep  = 10; // GeV
   // mHstart = 1; // GeV
   
   // for (int mHi = 0; mHi < 100; mHi++)
   //   {
   //     OSinput DS2l(4.40, 80.385, 91.1876, 125.6, 173.5);
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

