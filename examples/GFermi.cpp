#include <boost/math/tools/roots.hpp>
#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

struct DiffGF
{

  OSinput oi;
  WW<OS>* dW;
  ZZ<OS>* dZ;
  size_t order;

  DiffGF(OSinput in_, size_t order_) : oi(in_), order(order_)
  {
    dW = new WW<OS>(oi, oi.MMZ());
    dZ = new ZZ<OS>(oi, oi.MMZ());
  }

  long double operator()(long double alpha)
  {
    // long double Gf0 = 1.16637e-5;
    long double Gf0 = pdg2014::Gf;
    // long double alphaS = 0.1184;
    long double alphaS = pdg2014::asMZ;
      
    long double dMyW = 1;
    long double dMyZ = 1;
    long double Gf[4];
    // Tree level
    Gf[0] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);

    // 1-loop level
    dMyW += alpha/4./Pi*dW->y10();
    dMyZ += alpha/4./Pi*dZ->y10();
    
    Gf[1] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
    // 2-loop level
    dMyW += alpha/4./Pi*alphaS/4./Pi*dW->y11();
    dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ->y11();
    
    Gf[2] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
    
    dMyW += pow(alpha/4./Pi,2)*dW->y20();
    dMyZ += pow(alpha/4./Pi,2)*dZ->y20();
    
    Gf[3] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
    return (Gf0 - Gf[order])*pow(10,2);
  }

};

class tolerance {
public:
  tolerance(long double eps) :
    _eps(eps) {
  }
  bool operator()(long double a, long double b) {
    return (fabs(b - a) <= _eps);
  }
private:
  long double _eps;
};


int main (int argc, char *argv[])
{
  try
    {
      // Disable TSIL warnings output
      fclose(stderr);

      // Compare with:
      // OSinput KVPhys(4.9, 80.385, 91.1876, 125.6, 173.5);
      OSinput KVPhys(4.9, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt);


      DiffGF dGF(KVPhys,0);
      tolerance tol = 0.000000000001;
      
      std::cout << dGF(0.) * dGF(1.) << std::endl; 

      std::pair<long double, long double> found = boost::math::tools::bisect(dGF, 1./137., 1./126., tol);
      std::cout << std::setprecision(10);
      std::cout << "==> x = [" << 1./found.first << ',' << 1./found.second << "]\n";

      // Now separate contributions

      DiffGF dGF00(KVPhys,0);
      DiffGF dGF10(KVPhys,1);
      DiffGF dGF11(KVPhys,2);
      DiffGF dGF20(KVPhys,3);

      std::pair<long double, long double> found00 = boost::math::tools::bisect(dGF00, 1./137., 1./126., tol);
      std::pair<long double, long double> found10 = boost::math::tools::bisect(dGF10, 1./137., 1./126., tol);
      std::pair<long double, long double> found11 = boost::math::tools::bisect(dGF11, 1./137., 1./126., tol);
      std::pair<long double, long double> found20 = boost::math::tools::bisect(dGF20, 1./137., 1./126., tol);

      std::cout << " 1/alpha = " << 1./found00.first << std::endl;
      std::cout << " (1,0)   + " << 1./found10.first - 1./found00.first << std::endl;
      std::cout << " (1,1)   + " << 1./found11.first - 1./found10.first << std::endl;
      std::cout << " (2,0)   + " << 1./found20.first - 1./found11.first << std::endl;
      std::cout << " 2-l EW  = " << 1./found20.first  << std::endl;

      // Now analytic solution

      alphaGF aGF  = alphaGF(KVPhys, KVPhys.MMZ());
      
      long double Gf = pdg2014::Gf;
      long double alphaSMZ = pdg2014::asMZ;
      std::complex<long double> dMyW,dMyZ;
      dMyW = 1;
      dMyZ = 1;


      // Tree level
      long double alF = real(sqrt(2)*Gf*KVPhys.MMW()/Pi*(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ)*dMyW);

      // WW dW    = WW(KVPhys, KVPhys.MMZ());
      // ZZ dZ    = ZZ(KVPhys, KVPhys.MMZ());


      long double ali[3];
      ali[0] = alF;
      ali[1] = alF*(1 + alF/4./Pi*aGF.a10());
      ali[2] = alF*(1 + alF/4./Pi*aGF.a10() + alF/4./Pi*alphaSMZ/4./Pi*aGF.a11());
      ali[3] = alF*(1 + alF/4./Pi*aGF.a10() + alF/4./Pi*alphaSMZ/4./Pi*aGF.a11() + pow(alF/4./Pi,2)*aGF.a20());
      
      std::cout << std::setprecision(8);
      std::cout << "alpha=alpha-Born" << std::endl;
      std::cout << "1/alpha = " << 1./ali[0]  << " : born" << std::endl      
                << "          " << 1./ali[1] - 1./ali[0] << " : EW " << std::endl
                << "          " << 1./ali[2] - 1./ali[1] << " : EW * QCD " << std::endl
                << "          " << 1./ali[3] - 1./ali[2] << " : EW * EW " << std::endl
                << "          " << std::endl
                << "          " << 1./ali[3] << " : Total " << std::endl;


      return 0;

      long double alphaTree = 1./137.234;

      long double alphaS = 0.1184;
      // \mu = Mt
      long double alphaMt  = 0.00779305;

      
      // long double alphaMZ = 1./127.944;
      // long double alphaMZ = 1./127.773;
      long double alphaMZ = 1./126.654;



      // 
      // Test relation betwee alpha and G Fermi
      // 
      // 
      
      
      WW<OS> dW    = WW<OS>(KVPhys, KVPhys.MMZ());
      ZZ<OS> dZ    = ZZ<OS>(KVPhys, KVPhys.MMZ());
      std::ofstream of("GF.dat");
      for (double alpha = 1./137.; alpha < 1./125.; alpha += 0.00001)
        {
          long double dMyW = 1;
          long double dMyZ = 1;
          long double Gf[4];
          // Tree level
          Gf[0] = alpha*Pi/sqrt(2)/KVPhys.MMW()/dMyW/(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ);
          
          // 1-loop level
          dMyW += alpha/4./Pi*dW.y10();
          dMyZ += alpha/4./Pi*dZ.y10();
          
          Gf[1] = alpha*Pi/sqrt(2)/KVPhys.MMW()/dMyW/(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ);
          
          
          // 2-loop level
          dMyW += alpha/4./Pi*alphaS/4./Pi*dW.y11();
          dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ.y11();
          
          Gf[2] = alpha*Pi/sqrt(2)/KVPhys.MMW()/dMyW/(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ);
          
          
          dMyW += pow(alpha/4./Pi,2)*dW.y20();
          dMyZ += pow(alpha/4./Pi,2)*dZ.y20();
          
          Gf[3] = alpha*Pi/sqrt(2)/KVPhys.MMW()/dMyW/(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ);

          long double Gf0 = 1.16637e-5;

          of << 1./alpha << " " 
             << Gf0 - Gf[0] << " " 
             << Gf0 - Gf[1] << " " 
             << Gf0 - Gf[2] << " " 
             << Gf0 - Gf[3] << std::endl;
         

          std::cout << 1./alpha << " " 
                    << Gf[0] << " " 
                    << Gf[1] << " " 
                    << Gf[2] << " " 
                    << Gf[3] << std::endl;
        }
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

