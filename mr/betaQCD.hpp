#ifndef __BETAQCD_HPP_
#define __BETAQCD_HPP_
#include <boost/numeric/odeint.hpp>
#include "constants.hpp"

using namespace boost::numeric::odeint;

typedef std::vector<long double> state_type;

class BetaQCD
{

  size_t loops;
  double nf;

public:
  BetaQCD( size_t loops_, double nf_) : loops(loops_), nf(nf_) { }
  
  double b0()
  {
    return loops > 0 ? 
      (11.-2./3.*nf)/4. : 0;
  }
  double b1()
  {
    return loops > 1 ? 
      (102.-38./3.*nf)/16. : 0;
  }
  double b2()
  {
    return loops > 2 ? 
      (2857./2.-5033./18.*nf + 325./54.*nf*nf)/64. : 0;
  }
  double b3()
  {
    return loops > 3 ? 
      (149753./6. + 3564.*Zeta3
      - (1078361./162. + 6508./27.*Zeta3)*nf
      + (50065./162. + 6472./81.*Zeta3)*nf*nf
       + 1093./729.*nf*nf*nf)/256. : 0;
  }
  void operator() (const state_type &a, state_type &dadt, const double t)
  {
    dadt[0] = -b0()*pow(a[0],2) - b1()*pow(a[0],3) - b2()*pow(a[0],4) - b3()*pow(a[0],5);
  }
};


class AlphaS
{
  double   mu0;
  size_t loops;
  size_t    nf;
public:
  AlphaS(size_t loops_ = 4, size_t nf_ = 5) : loops(loops_), nf(nf_)
  { }
  double operator()(long double mu2)
  {
    state_type a0(1);

    // Starting value
    a0[0] = 0.1184/Pi; 
            
    BetaQCD beta4l5nf(loops,nf);
    
    double mu0 = pow(91.1876,2);
    
    double lEnd = log(mu2/mu0);
    
    controlled_stepper_standard< stepper_rk5_ck< state_type > >
      controlled_rk5( 1E-6 , 1E-7 , 1.0 , 1.0 );
 
    integrate_adaptive( controlled_rk5 ,
                        beta4l5nf, a0, 0.0, lEnd, 0.01  );
    
    return a0[0]*Pi; 

  }
};

#endif  // __BETAQCD_HPP_
