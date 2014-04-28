#ifndef __BETAQCD_HPP_
#define __BETAQCD_HPP_
#include <boost/numeric/odeint.hpp>
// #include "LambertW.h"
using namespace boost::numeric::odeint;
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

// class BetaQCD
// {

//   size_t loops;
//   double nf;

// public:
//   BetaQCD( size_t loops_, double nf_) : loops(loops_), nf(nf_) { }
  
//   double b0()
//   {
//     return loops > 0 ? 
//       (11.-2./3.*nf)/4. : 0;
//   }
//   double b1()
//   {
//     return loops > 1 ? 
//       (102.-38./3.*nf)/16. : 0;
//   }
//   double b2()
//   {
//     return loops > 2 ? 
//       (2857./2.-5033./18.*nf + 325./54.*nf*nf)/64. : 0;
//   }
//   double b3()
//   {
//     return loops > 3 ? 
//       (149753./6. + 3564.*Zeta3
//       - (1078361./162. + 6508./27.*Zeta3)*nf
//       + (50065./162. + 6472./81.*Zeta3)*nf*nf
//        + 1093./729.*nf*nf*nf)/256. : 0;
//   }
//   void operator() (const state_type &a, state_type &dadt, const double t)
//   {
//     dadt[0] = -b0()*pow(a[0],2) - b1()*pow(a[0],3) - b2()*pow(a[0],4) - b3()*pow(a[0],5);
//   }
// };

class AlphaS
{
  double   mu0;
  size_t loops;
public:
  AlphaS(size_t loops_) : loops(loops_)
  { }
  double operator()(double mu2)
  {
    state_type a0(1);
    a0[0] = 0.1184/Pi; // start at x=1.0, p=0.0
    
    // Starting value
    
    std::vector<double> times;
    std::vector<state_type> x_t_vec;
    
    BetaQCD beta4l5nf(loops,5);
    
    double mu0 = pow(91.1876,2);
    // double mu2  = pow(mu,2);
    
    double lEnd = log(mu2/mu0);
    
    size_t steps = integrate(  beta4l5nf, 
                               a0 , 0.0 , lEnd , 
                               back_inserter( times ) ,
                               back_inserter( x_t_vec ) );
    return x_t_vec[steps][0]*Pi;
    
  }
};

// Exact two-loop solution with the 
// help of Lambert W-function
// class AlphaS2l
// {
// public:
//   double operator()(double mu2)
//   {
//     double nf = 5.;
//     double be0 = 11. - 2./3.*nf;
//     double be1 = 102. - 38./3.*nf; 
//     double B1  = be1/be0/be0;
//     double mmu0 = pow(91.1876,2);

//     // alphaS at M_Z
//     double asMZ = 0.1184;
//     double Lam2=mmu0*exp(-4*Pi/be0/asMZ);
//     Lam2=pow(0.213,2);

//     std::cout << "Lambda QCD (MeV)= " << sqrt(Lam2)*1000. << std::endl;
//     return -4*Pi/be0/B1/(1 + LambertW(-1, -exp(-(1+log(mu2/Lam2)/B1))));
//   }
// };

#endif  // __BETAQCD_HPP_
