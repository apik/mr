#include <boost/math/tools/roots.hpp>
#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

#include <boost/numeric/interval.hpp>



struct DiffLambda
{

  // OSinput oi;
  // WW<OS>* dW;
  // ZZ<OS>* dZ;
  // size_t order;
  RunUpto ru;

  DiffLambda(RunUpto ru_) : ru(ru_)
  {
    // dW = new WW<OS>(oi, oi.MMZ());
    // dZ = new ZZ<OS>(oi, oi.MMZ());
    
  }

  long double operator()(long double mu)
  {
    std::cout << " lam = " << ru.lambda(mu) << std::endl;
    return ru.lambda(mu);
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




long double alpha(long double mu)
{
  long double MZ  = 91.1876;
  long double aMZ = 1./127.944;
  return aMZ/(1+11./6./Pi*aMZ*log(mu/MZ));
}

int main (int argc, char *argv[])
{
  try
    {

      // Disable TSIL warnings output
      fclose(stderr);
      long double MMt,MMW,MMZ,MMH,alphaMt,alphaS,alphaSMt;

      // Scale inv test
      long double a = 1.;
      // Compare with:
      alphaMt   = 1./128.175;
      alphaSMt   = 0.1079;
      
      std::vector<OSinput> sv;

      sv.push_back(OSinput(0, 80.384, 91.1876, 125.66, 173.10)); // A
      // sv.push_back(OSinput(a*4.4, a*80.399, a*91.1876, a*126, a*173.9)); // B
      // sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125, a*172.9)); // C

      // sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125, a*274.1)); // C
      // sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125, a*374.1)); // C
      // sv.push_back(OSinput(a*0, a*80.399, a*91.1876, a*125, a*474.1)); // C

      sv.push_back(OSinput(0, pdg2014::MW, pdg2014::MZ, pdg2014::MH, pdg2014::Mt)); // 

      RunUpto rr(sv[0], pdg2014::aMZ, 0.1184);

      std::cout << " LAMMMMMMMMMMMM " << rr.lambda(10000) << std::endl;
      
      DiffLambda dl(rr);

      tolerance tol = 0.1;
      
      std::pair<long double, long double> found = boost::math::tools::bisect(dl, 10e5, 10e25, tol);
      std::cout << std::setprecision(10);
      std::cout << "==> x = [" << found.first << ',' << found.second << "]\n";

      
      std::cout << "\n\nPDG 2014" << std::endl; 
      RunUpto rrPDG2014(sv[1]);
      

      // Degrassi table 3

      OSinput oi(0, 80.384, 91.1876, 125.66, 173.10);
      
      AlphaS       as;
      AlphaQEDQCD  al;
      long double ms = oi.Mt();
      
      WW<OS> w(oi, pow(ms,2));
      ZZ<OS> z(oi, pow(ms,2));
      HH<OS> h(oi, pow(ms,2));
      tt<OS> t(oi, pow(ms,2));
      
      // long double aQCD = as(ms)/4./Pi;
      // long double aEW  = al.QED(ms)/4./Pi;

      long double aQCD = as(ms)/4./Pi;
      long double aEW  = al.QED(ms)/4./Pi;

      long double gsD  = pow(0.65294,2);
      long double gpsD = pow(0.34972,2);
      // aEW = gsD*gpsD/(gsD+gpsD)/16./Pi/Pi;
      
      std::cout << "1/al = " << 1./(aEW*4.*Pi) << std::endl;
      std::cout << "aS = " << (aQCD*4.*Pi) << std::endl;

      
      // 2-loop EW corrections
      long double dWplus1 = 1 + aEW*w.y10() + 0*(aEW*aQCD*w.y11() + aEW*aEW*w.y20());
      long double dZplus1 = 1 + aEW*z.y10() + 0*(aEW*aQCD*z.y11() + aEW*aEW*z.y20());
      long double dHplus1 = 1 + aEW*h.y10() + 0*(aEW*aQCD*h.y11() + aEW*aEW*h.y20());
      long double dtplus1 = 1 + aEW*t.y10() + 0*(aEW*aQCD*t.y11() + aEW*aEW*t.y20())
        + aQCD*t.y01()+ 0*(aQCD*aQCD*t.y02()+ aQCD*aQCD*aQCD*t.y03());
      
      long double Gf = pdg2014::Gf;
      
      long double gg = pow(2.,5./2.)*Gf*oi.MMW()*dWplus1;
      long double gg_ggp = pow(2.,5./2.)*Gf*oi.MMZ()*dZplus1;
      
      long double a1 = 5./3.*(gg_ggp - gg)/16./Pi/Pi;
      long double a2 = gg/16./Pi/Pi;
      long double aS = aQCD;
      long double ayt = pow(2.,3./2.)*Gf*oi.MMt()*pow(dtplus1,2)/16./Pi/Pi;
      long double alam = Gf/sqrt(2.)*oi.MMH()*dHplus1/16./Pi/Pi;
      
      
      std::cout << " At matching scale mu = " << ms << std::endl;
      std::cout << " g1 = " << sqrt(3./5.*a1)*4*Pi << std::endl;
      std::cout << " g2 = " << sqrt(a2)*4*Pi << std::endl;
      std::cout << " g3 = " << sqrt(aS)*4*Pi << std::endl;
      std::cout << " yt = " << sqrt(ayt)*4*Pi << std::endl;
      std::cout << " lam = " << alam*16*Pi*Pi << std::endl;
      
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

