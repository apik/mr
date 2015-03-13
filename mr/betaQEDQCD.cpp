#include "betaQEDQCD.hpp"
#include "boost/numeric/odeint/integrate/integrate_adaptive.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/numeric/odeint/stepper/controlled_runge_kutta.hpp"

using namespace boost::numeric::odeint;
typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

//         QCD  &  QED
std::pair<double,double> runQEDQCD(long double aStart, long double asStart, long double muStart, long double muEnd, size_t ng)
{
  state_type alpha4pi(2);
  
  // Starting value
  alpha4pi[0] = asStart/4./Pi;
  alpha4pi[1] =  aStart/4./Pi; 
  
  double lEnd = 2.*log(muEnd/muStart);

  
  // std::cout << "End time: " << fabs(lEnd) << std::endl;
  BetaQEDQCD bEMQCD(ng, lEnd < 0);
    
  double abs_err = 1E-6;
  double rel_err = 1E-7;
  double a_x = 1.0;
  double a_dxdt = 1.0;
  
  // OdeInt_v.1
  // controlled_stepper_standard< stepper_rk5_ck< state_type > >
  //   controlled_rk5( 1E-6 , 1E-7 , 1.0 , 1.0 );
  
  // integrate_adaptive( controlled_rk5 ,
  //                     bEMQCD, alpha4pi, 0.0, fabs(lEnd), 0.01  );
  
  controlled_stepper_type 
    controlled_stepper(default_error_checker< double , range_algebra , default_operations >
                       ( abs_err , rel_err , a_x , a_dxdt ) );
  
  integrate_adaptive( controlled_stepper ,
                      bEMQCD, alpha4pi, 0.0, fabs(lEnd), 0.01  );
  
  return std::make_pair(4.*Pi*alpha4pi[0],4.*Pi*alpha4pi[1]);            // We return a_S
}



