#include <iostream>
#include <cmath>
#include "mr.hpp"
#include "mfitter.hpp"

#include "Minuit/FunctionMinimum.h"
#include "Minuit/MnUserParameterState.h"
#include "Minuit/MnPrint.h"
#include "Minuit/MnMigrad.h"
#include "Minuit/MnMinos.h"
#include "Minuit/MnContours.h"
#include "Minuit/MnPlot.h"
#include "Minuit/MinosError.h"
#include "Minuit/ContoursError.h"

int main (int argc, char *argv[])
{
  try
    {

      std::vector<long double> ob = observables(0.35761,
                                                0.64822,
                                                1.1666,
                                                0.93558,
                                                0.12711,
                                                132.03,
                                                pow(173.10,2)
                                                );


      std::cout << "MW = " << ob[0] << std::endl
                << "MZ = " << ob[1] << std::endl
                << "MH = " << ob[2] << std::endl
                << "Mt = " << ob[3] << std::endl
                << "Gf = " << ob[4] << std::endl;


      std::cout << "\n\n Clear input: " << std::endl;
      MSinput mmx = MSinput::fromConsts(pow(173.10,2), 132.03, 0.12711, 0, 0.93558, 0.64822, 0.35761);


      std::cout << mmx.mW() << std::endl;
      std::cout << mmx.mZ() << std::endl;
      std::cout << mmx.mH() << std::endl;
      std::cout << mmx.mt() << std::endl;
      std::cout << mmx.vev() << std::endl;


      // std::vector<double> pos(6);

      std::vector<double> meas(6);
      std::vector<double> var(6);

      // MW
      meas[0] = 80.6092;
      var[0]  =  pow(0.015,2);
      
      // MZ
      meas[1] = 91.3923;
      var[1]  =  pow(0.0021,2);
      
      // MH
      meas[2] = 125.912;
      var[2]  =   pow(0.4,2);
      
      // Mt
      meas[3] = 168.913;
      var[3]  = pow(0.51,2);
      
      // Gf
      meas[4] = 1.15658e-5;
      var[4]  = pow(0.0000006,2);

      // alphaS
      meas[5] = 0.108301;
      var[5]  = pow(0.0006,2);

  
      // // MW
      // meas[0] = 80.385;
      // var[0]  =  pow(0.015,2);
      
      // // MZ
      // meas[1] = 91.1876;
      // var[1]  =  pow(0.0021,2);
      
      // // MH
      // meas[2] = 125.7;
      // var[2]  =   pow(0.4,2);
      
      // // Mt
      // meas[3] = 173.21;
      // var[3]  = pow(0.51,2);
      
      // // Gf
      // meas[4] = 1.1663787;
      // var[4]  = pow(0.0000006,2);

      // // alphaS
      // meas[5] = 0.1185;
      // var[5]  = pow(0.0006,2);

      // create FCN function  
      MfitterFcn theFCN(meas, var, 173.4);

      // // starting values for parameters
      // std::vector<double> init_par; 
      // init_par.push_back(0.35761);
      // init_par.push_back(0.64822);
      // init_par.push_back(1.1666);
      // init_par.push_back(0.93558); 
      // init_par.push_back(0.12711); 
      // init_par.push_back(132.03); 

      // std::cout << "CHI= " << theFCN(init_par) << std::endl;

      // return 0;
      
      long double gpIN  = 0.35761;
      long double gIN   = 0.64822;
      long double gsIN  = 1.1666;
      long double ytIN  = 0.93558;
      long double lamIN = 0.12711;
      long double mu0IN = 132.03;
                  
      MnUserParameters upar;
      upar.add("gp",  gpIN, 0.1);
      upar.add("g",   gIN,  0.1);
      upar.add("gs",  gsIN, 0.1);
      upar.add("yt",  ytIN, 0.1);
      upar.add("lam", lamIN,0.1);
      upar.add("mu0", mu0IN,0.1);

      
      
      // access parameter by name to set limits...
      upar.setLimits("gp", gpIN-0.2, gpIN+0.2);
      upar.setLimits("g", gIN-0.2, gIN+0.2);
      upar.setLimits("gs", gsIN-0.3, gsIN+0.3);
      upar.setLimits("yt", ytIN-0.2, ytIN+0.2);
      upar.setLimits("lam", lamIN-0.5, lamIN+0.5);
      upar.setLimits("mu0", mu0IN-20, mu0IN+20);

      MnMigrad migrad(theFCN, upar);
      
      // // fix a parameter...
      // migrad.fix("gp");
      // migrad.fix("g");
      // migrad.fix("gs");
      // migrad.fix("yt");
      // migrad.fix("mu0");      

      // // ... and minimize
      FunctionMinimum min = migrad();
      
      // output
      std::cout<<"\n\n\n minimum: "<<min<<std::endl;


      // create MINOS error factory
      // MnMinos minos(theFCN, min);
      
      // {
      //   // 1-sigma MINOS errors (minimal interface)
      //   std::pair<double,double> e0 = minos(0);
      //   std::pair<double,double> e1 = minos(1);
      //   std::pair<double,double> e2 = minos(2);
      //   std::pair<double,double> e3 = minos(3);
      //   std::pair<double,double> e4 = minos(4);
      //   std::pair<double,double> e5 = minos(5);

      //   // output
      //   std::cout<<"1-sigma minos errors: "<<std::endl;
      //   std::cout<<"par0: "<<min.userState().value("gp")<<" "<<e0.first<<" "<<e0.second<<std::endl;
      //   std::cout<<"par1: "<<min.userState().value("g")<<" "<<e1.first<<" "<<e1.second<<std::endl;
      //   std::cout<<"par2: "<<min.userState().value("gs")<<" "<<e2.first<<" "<<e2.second<<std::endl;
      //   std::cout<<"par3: "<<min.userState().value("yt")<<" "<<e3.first<<" "<<e3.second<<std::endl;
      //   std::cout<<"par4: "<<min.userState().value("lam")<<" "<<e4.first<<" "<<e4.second<<std::endl;
      //   std::cout<<"par5: "<<min.userState().value("mu0")<<" "<<e5.first<<" "<<e5.second<<std::endl;
        
      // }

      
      // // starting values for initial uncertainties
      // std::vector<double> init_err; 
      // init_err.push_back(0.01); 
      // init_err.push_back(0.01); 
      // init_err.push_back(0.01);
      // init_err.push_back(0.01); 
      // init_err.push_back(0.01); 
      // init_err.push_back(0.1);
      
      // // create minimizer (default constructor)
      // VariableMetricMinimizer theMinimizer;
      
      // // minimize
      // FunctionMinimum min = theMinimizer.minimize(theFCN, init_par, init_err);
      
    // output
    // std::cout<<"minimum: "<<min<<std::endl;


    }
  catch (std::exception &p) 
    {
      std::cerr << p.what() << std::endl;
      return 1;
    }
  
  return 0;
}

