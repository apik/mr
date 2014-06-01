#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"
#include "CRunDec.h"


// CRunDec  input
#define asMz 0.1184
#define Mz 91.1876
#define Mt 173.5
#define Mb 4.4
#define Mc 1.5
#define muc 1.279
#define mub 4.163
#define Mtau 1.777

long double as(long double mu2)
{
  long double mu = sqrt(mu2);
  CRunDec crundec;
  crundec.nfMmu[0].nf = 5;
  crundec.nfMmu[0].Mth = Mb;
  crundec.nfMmu[0].muth = Mb;
  crundec.nfMmu[1].nf = 6;
  crundec.nfMmu[1].Mth = Mt;
  crundec.nfMmu[1].muth = Mt;
  
  if (mu < Mz) 
    return crundec.AlH2AlL(asMz, Mz,crundec.nfMmu,mu,4);
  else
    return crundec.AlL2AlH(asMz, Mz,crundec.nfMmu,mu,4);
}

long double a(long double mm)
{

  long double alphaMZ = 1./127.916;
  long double mme1,mme2,mme3,mmW,mmu3;

  mme1 = pow(0.51/1000.,2);
  mme2 =pow(105.6583/1000.,2);
  mme3 =pow(1777.03/1000.,2);
  mmu3 =pow(173.5,2);
  mmW = pow(80.385,2);
  long double mZpole = 91.1876;
  long double HADRONS = 0.027690;
  // return alpha0*(7.*log(mmW/mm)-2./3. - 4./3.*( log(mme1/mm) + log(mme2/mm) + log(mme3/mm) )
  //                        - 16/9*log(mmu3/mm)
  //                        - 44/9*(log(mmW/mm) - 5/3)
  //                        )
  //   + HADRONS ;

  long double b1=-22./3+4./3.*3.+1./6;
  long double b1p=(4./3*3.+1./10)*5./3.;
  return alphaMZ/(1-(b1+b1p)/(2*Pi) * alphaMZ * log(sqrt(mm)/mZpole));


}

int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS,alphaSMt;

      // Compare with:
      SMinput KVPhys(4.40, 80.385, 91.1876, 125.66, 173.5);
      


      long double alphaSMZ = 0.1184;
      // \mu = Mt
      alphaMt  = 0.00779305;
      alphaSMt   = 0.1184;

      
      // \mu = Mb
      long double alphaMb  = 0.00784257;
      long double alphaSMb   = 0.1905;

      

      AlphaS as;
      
      long double alphaMZ = 1./137.035999;
      
      
   //    // Yukawa top
   //    tt dMt  = tt(KVPhys, KVPhys.MMt());
   //    std::cout << "[ Top quark ]" << std::endl;
   //    std::cout << "Mh= " << KVPhys.MH()  << std::endl;
   //    std::cout << "as(MMt) = " << as(KVPhys.MMt()) << std::endl;          
   //    std::cout << "\t1-loop \\alpha         " << alphaMt/4./Pi*dMt.my10() << std::endl;
   //    std::cout << "\t1-loop \\alpha_S       " << alphaSMt/4./Pi*dMt.my01() << std::endl;
   //    std::cout << "\t2-loop \\alpha*\\alpha_S" << alphaMt/4./Pi*alphaSMt/4./Pi*dMt.my11() << std::endl;
   //    std::cout << "\t2-loop \\alpha^2       " << pow(alphaMt/4./Pi,2)*dMt.my20() << std::endl;


   //    // Yukawa bottom
   //    bb dMb  = bb(KVPhys, KVPhys.MMb());
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


      // Plot Yukawa top
      Plot3 plotYt("yt", "Yukawa Top", "\\mu", "Yukawa-top", "1-loop,a^1", "2-loop,a*a_S", "2-loop,a^2");
      long double mHstep  = 50; // GeV
      long double mHstart = 100; // GeV
      
   for (int mHi = 0; mHi < 19; mHi++)
     {
       SMinput DS2l(4.40, 80.385, 91.1876, 125.66, 173.5);
       long double mu = mHstart + mHi*mHstep;
       tt dtY  = tt(DS2l, mu*mu);          
       
       std::cout << "mu= " << mu << ",  as = " << as(mu) << ",  dyt = " << dtY.my11().real() << std::endl;
       plotYt.add(mu,a(mu)/4./Pi*dtY.my10().real(),
                  // a(mu)/4./Pi*as(mu)/4./Pi*
                  dtY.my11().real(),
                  // pow(a(mu)/4./Pi,2)*
                  dtY.my20().real());

     }


   // Plot Yukawa bottom
   Plot3 plotYb("yb", "Yukawa Bottom", "\\mu", "Yukawa bottom", "a^1", "a*as","a^2");
   mHstep  = 10; // GeV
   mHstart = 1; // GeV
   
   for (int mHi = 0; mHi < 100; mHi++)
     {
       SMinput DS2l(4.40, 80.385, 91.1876, 125.6, 173.5);
       long double mu = mHstart + mHi*mHstep;
       bb dbY  = bb(DS2l, mu*mu);          
       
       plotYb.add(mu,a(mu)/4./Pi*dbY.my10().real(),
                  a(mu)/4./Pi*as(mu)/4./Pi*dbY.my11().real(),
                  pow(a(mu)/4./Pi,2)*dbY.my20().real());
     }

    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

