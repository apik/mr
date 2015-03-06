#include <boost/math/tools/roots.hpp>
#include <iostream>
#include "mr.hpp"
#include "gnuplot.hpp"

// #include <boost/numeric/interval.hpp>

#include "errint.hpp"


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




// long double alpha(long double mu)
// {
//   long double MZ  = 91.1876;
//   long double aMZ = 1./127.944;
//   return aMZ/(1+11./6./Pi*aMZ*log(mu/MZ));
// }

int main (int argc, char *argv[])
{
  try
    {

      // Disable TSIL warnings output
      fclose(stderr);

      OSinputErr oie(4.9, we(pdg2014::MW,0.555555555), pdg2014::MZ, pdg2014::MH, pdg2014::Mt);

      std::cout  << "MMW = ( " << oie.MMW().lower() << ", " << oie.MMW().upper() << " )\n"
                << "Empty: " << boost::numeric::width(oie.MW()) << std::endl;
      
      WWerr wwe(oie,1);
      
      
    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

