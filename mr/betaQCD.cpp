#include "betaQCD.hpp"
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>

using namespace boost::numeric::odeint;
typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;


double run(long double asStart, long double muStart, long double muEnd, size_t NF, size_t loops)
{
  state_type as4pi(1);
  
  // Starting value
  as4pi[0] = asStart/4./Pi; 
  
  double lEnd = 2.*log(muEnd/muStart);

  
  std::cout << "End time: " << fabs(lEnd) << std::endl;
  BetaQCD beta4l5nf(loops, NF, lEnd < 0);
    
  double abs_err = 1E-6;
  double rel_err = 1E-7;
  double a_x = 1.0;
  double a_dxdt = 1.0;
  
  // OdeInt_v.1
  // controlled_stepper_standard< stepper_rk5_ck< state_type > >
  //   controlled_rk5( 1E-6 , 1E-7 , 1.0 , 1.0 );
  
  // integrate_adaptive( controlled_rk5 ,
  //                     beta4l5nf, as4pi, 0.0, fabs(lEnd), 0.01  );

  // OdeInt_v.2
  controlled_stepper_type 
    controlled_stepper(default_error_checker< double , range_algebra , default_operations >
                       ( abs_err , rel_err , a_x , a_dxdt ) );
  integrate_adaptive( controlled_stepper , beta4l5nf , as4pi , 0.0 , fabs(lEnd) , 0.01 );

  
  return as4pi[0]*4.*Pi;            // We return a_S
}



// Based on expression eq.63 from
// Chetyrkin, Kuhn, Sturm 
// Nucl.Phys. B744 (2006) 121-135
// arXiv:hep-ph/0512060
long double as5nf2as6nf(long double M, long double mu, long double as, size_t nl, size_t ord)
{
  
  long double as4Pi = as/4./Pi;
  
  long double Lmu2MM = 2.*log(mu/M);
  
  // z[as(nf=5)] = as(nf=6)/as(nf=5)
  long double z = 1;

  if(ord > 0)
    
    z += as4Pi*Lmu2MM * ( 2./3. );
  
  if(ord > 1)
    {
      
      z +=  + pow(as4Pi,2) * ( 14./3. );
      
      z +=  + pow(as4Pi,2)*Lmu2MM * ( 38./3. );
      
      z +=  + pow(as4Pi,2)*pow(Lmu2MM,2) * ( 4./9. );
      
    }
  
  if(ord > 2)
    {
      
      z +=  + pow(as4Pi,3) * ( 58933./1944. + 80507./432.*Zeta3 + 128./3.*Zeta2 + 128./9.*log(2)*Zeta2 );
      
      z +=  + pow(as4Pi,3)*nl * (  - 2479./486. - 64./9.*Zeta2 );
      
      z +=  + pow(as4Pi,3)*Lmu2MM * ( 8941./27. );
      
      z +=  + pow(as4Pi,3)*Lmu2MM*nl * (  - 409./27. );
      
      z +=  + pow(as4Pi,3)*pow(Lmu2MM,2) * ( 511./9. );

      z +=  + pow(as4Pi,3)*pow(Lmu2MM,3) * ( 8./27. );
    }

  if(ord > 3)
    {

      z +=  + pow(as4Pi,4) * ( 592371712./382725. - 21814592./2835.*
                               a5 - 25433192./1701.*a4 + 40596749./5670.*Zeta5 + 71102219./8505.
                               *Zeta4 + 2408412383./340200.*Zeta3 + 11153936./1215.*Zeta2 - 46048./
                               27.*Zeta2*Zeta3 - 18636934./2835.*log(2)*Zeta4 - 131456./81.*log(2)*Zeta2 + 
                               5826074./1701.*pow(log(2),2)*Zeta2 - 5453648./8505.*pow(log(2),3)*Zeta2
                               - 3179149./5103.*pow(log(2),4) + 2726824./42525.*pow(log(2),5) );
      
      z +=  + pow(as4Pi,4)*nl * (  - 1773073./2916. - 692./81.*a4 - 
                                   460./9.*Zeta5 + 697709./648.*Zeta4 - 4756441./3888.*Zeta3 - 71296./
                                   81.*Zeta2 - 5632./81.*log(2)*Zeta2 + 1709./81.*pow(log(2),2)*Zeta2 - 173./
                                   486.*pow(log(2),4) );
      
      z +=  + pow(as4Pi,4)*pow(nl,2) * ( 140825./5832. + 76./27.*Zeta3
                                         + 1664./81.*Zeta2 );
      
      z +=  + pow(as4Pi,4)*Lmu2MM * ( 21084715./2916. + 2922161./648.*Zeta3
                                 + 8960./9.*Zeta2 + 8960./27.*log(2)*Zeta2 );
      
      z +=  + pow(as4Pi,4)*Lmu2MM*nl * (  - 1140191./1458. - 132283./324.*
                                     Zeta3 - 6016./27.*Zeta2 - 512./27.*log(2)*Zeta2 );
      
      z +=  + pow(as4Pi,4)*Lmu2MM*pow(nl,2) * ( 1679./729. + 256./27.*Zeta2);
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,2) * ( 94078./27. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,2)*nl * (  - 18230./81. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,2)*pow(nl,2) * ( 493./81. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,3) * ( 28298./81. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,3)*nl * (  - 428./27. );
      
      z +=  + pow(as4Pi,4)*pow(Lmu2MM,4) * ( 16./81. );
      
    }

  if(ord > 4)
    std::cerr << "WARNING: Only 4-loop decoupling relations available."<< std::endl;

  return as*z;
}
