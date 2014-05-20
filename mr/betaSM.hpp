#ifndef __BETASM_HPP__
#define __BETASM_HPP__

#include <iomanip>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include "constants.hpp"

using namespace boost::numeric::odeint;

typedef std::vector<long double> state_type;

typedef std::vector<size_t> index_t;

struct index_cmp_t : std::binary_function<index_t, index_t, bool> {
  bool operator ()(index_t const& a, index_t const& b) const {
    for (index_t::size_type i = 0; i < a.size(); ++i)
      if (a[i] != b[i])
        return a[i] < b[i];
    return false;
  }
};



index_t CouplingsPowers(size_t, size_t, size_t, size_t, size_t, size_t, size_t);





class BetaSMFull
{

  int NG;
  // Maximal numer of loops available for beta-functions
  static const size_t maxBetaOrder = 3;

  static const size_t maxPoco = maxBetaOrder + 2;

  std::map<index_t, long double, index_cmp_t> be1;
  std::map<index_t, long double, index_cmp_t> be2;
  std::map<index_t, long double, index_cmp_t> be3;
  std::map<index_t, long double, index_cmp_t> be4;
  std::map<index_t, long double, index_cmp_t> be5;
  std::map<index_t, long double, index_cmp_t> be6;
  std::map<index_t, long double, index_cmp_t> be7;
  
  // int pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam;
  int pocoa1, pocoa2, pocoa3, pocoa4, pocoa5, pocoa6, pocoa7;
  size_t maxPower;

  void add(std::map<index_t, long double, index_cmp_t>&, size_t, size_t, size_t, size_t, size_t, size_t, size_t, size_t, long double);

public:

  // BetaSMFull(int pocoa1_, int pocoa2_, int pocoas_, int pocoat_, int pocoab_, int pocoatau_, int pocolam_);
  BetaSMFull(int, int, int, int, int, int, int, int);

  void operator() (const state_type &, state_type &, const double);
  long double betaQCD(long double);
};







// 
// Couplings class
// 
// template<size_t pocoGauge, size_t pocoYukawa, size_t pocoLam>
template < int pocoa1, int pocoa2, int pocoas, int pocoat, int pocoab, int pocoatau, int pocolam > 
class CouplingsSM
{
  
  double      mu0;
  state_type aSM0;
  size_t       NG;

public:
  CouplingsSM(double a1, double a2, double as, double at, double ab, double atau, double lam, double mu0_, size_t NG_) : mu0(mu0_), NG(NG_)
    {

      aSM0.push_back(a1);
      aSM0.push_back(a2);
      aSM0.push_back(as);
      aSM0.push_back(at);
      aSM0.push_back(ab);
      aSM0.push_back(atau);
      aSM0.push_back(lam);
      
      const size_t fw = 20;
      std::cout << std::scientific << std::setprecision(fw-10);
      std::cout << "# [Initial values:]" << std::endl;
      std::cout << "\t  a1 = " << std::setw(fw) << a1 
                << "   a2 = " << std::setw(fw) << a1 
                << "   as = " << std::setw(fw) << as <<std::endl; 
      std::cout << "\t  at = " << std::setw(fw) << at 
                << "   ab = " << std::setw(fw) << ab 
                << " atau = " << std::setw(fw) << atau <<std::endl; 
      std::cout << "\talam = " << std::setw(fw) << lam
                << "   MU = " << std::setw(fw) << sqrt(mu0)
                << "   NG = " << std::setw(fw) << NG <<std::endl; 

    }

  state_type operator()(long double mu2)
  {
        
    BetaSMFull be(pocoa1, pocoa2, pocoas, pocoat, pocoab, pocoatau, pocolam, NG);
    
    double lEnd = log(mu2/mu0);
    
    controlled_stepper_standard< stepper_rk5_ck< state_type > >
      controlled_rk5( 1E-6 , 1E-7 , 1.0 , 1.0 );
 
    state_type aSM(aSM0);
    integrate_adaptive( controlled_rk5 ,
                        be, aSM, 0.0, lEnd, 0.01
                        );
    return aSM; 
  }
};

#endif  // __BETASM_HPP__

