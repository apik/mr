#include <boost/math/tools/roots.hpp>
#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

struct DiffGF
{

  OSinput oi;
  WW<OS>* dW;
  ZZ<OS>* dZ;

  DiffGF(OSinput in_) : oi(in_)
  {
    dW = new WW<OS>(oi, oi.MMZ());
    dZ = new ZZ<OS>(oi, oi.MMZ());
  }

  long double operator()(long double alpha)
  {
    long double Gf0 = 1.16637e-5;
    long double alphaS = 0.1184;
      
    long double dMyW = 1;
    long double dMyZ = 1;
    long double Gf[4];
    // Tree level
    Gf[0] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
          
    // 1-loop level
    dMyW += alpha/4./Pi*dW->my10().real();
    dMyZ += alpha/4./Pi*dZ->my10().real();
          
    Gf[1] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
          
          
    // 2-loop level
    dMyW += alpha/4./Pi*alphaS/4./Pi*dW->my11().real();
    dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ->my11().real();
          
    Gf[2] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
          
          
    dMyW += pow(alpha/4./Pi,2)*dW->my20().real();
    dMyZ += pow(alpha/4./Pi,2)*dZ->my20().real();
          
    Gf[3] = alpha*Pi/sqrt(2)/oi.MMW()/dMyW/(1-oi.MMW()/oi.MMZ()*dMyW/dMyZ);
    
    return (Gf0 - Gf[3])*pow(10,2);
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
      OSinput KVPhys(4.40, 80.385, 91.1876, 125.6, 173.5);


      DiffGF dGF(KVPhys);
      tolerance tol = 0.0000000001;
      
      std::cout << dGF(0.) * dGF(1.) << std::endl; 

      std::pair<long double, long double> found = boost::math::tools::bisect(dGF, 1./137., 1./126., tol);
      std::cout << std::setprecision(10);
      std::cout << "==> x = [" << 1./found.first << ',' << 1./found.second << "]\n";
      
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
          dMyW += alpha/4./Pi*dW.my10().real();
          dMyZ += alpha/4./Pi*dZ.my10().real();
          
          Gf[1] = alpha*Pi/sqrt(2)/KVPhys.MMW()/dMyW/(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ);
          
          
          // 2-loop level
          dMyW += alpha/4./Pi*alphaS/4./Pi*dW.my11().real();
          dMyZ += alpha/4./Pi*alphaS/4./Pi*dZ.my11().real();
          
          Gf[2] = alpha*Pi/sqrt(2)/KVPhys.MMW()/dMyW/(1-KVPhys.MMW()/KVPhys.MMZ()*dMyW/dMyZ);
          
          
          dMyW += pow(alpha/4./Pi,2)*dW.my20().real();
          dMyZ += pow(alpha/4./Pi,2)*dZ.my20().real();
          
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

