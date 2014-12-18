// Example of YAML input

// MW:
//   mean: 80.385
//   error: 0.015
// MZ:
//   mean: 91.1876
//   error: 0.0021
// MH:
//   mean: 125.7
//   error:  0.4
// Mt:
//   mean: 173.21
//   error:  1.22                  # (0.51)(0.71)
// Gf:
//   mean:  0.000011663787
//   error: 0.000000000006
// as:
//   mean:  0.108057
//   error: 0.000496
// mu: 173.21




#include <iostream>
#include <fstream>
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

#include <yaml-cpp/yaml.h>


int main (int argc, char *argv[])
{
  try
    {

      std::vector<double> meas(6);
      std::vector<double> var(6);

      YAML::Node config;
      if(argc == 2)
        config = YAML::LoadFile(argv[1]);
      else
        config = YAML::LoadFile("fit.yaml");

      if (config["MW"]) {
        meas[0] = config["MW"]["mean"].as<double>();
        var[0] = pow(config["MW"]["error"].as<double>(),2);
      }
      if (config["MZ"]) {
        meas[1] = config["MZ"]["mean"].as<double>();
        var[1] = pow(config["MZ"]["error"].as<double>(),2);
      }
      if (config["MH"]) {
        meas[2] = config["MH"]["mean"].as<double>();
        var[2] = pow(config["MH"]["error"].as<double>(),2);
      }
      if (config["Mt"]) {
        meas[3] = config["Mt"]["mean"].as<double>();
        var[3]  = pow(config["Mt"]["error"].as<double>(),2);
      }
      if (config["Gf"]) {
        meas[4] = config["Gf"]["mean"].as<double>();
        var[4]  = pow(config["Gf"]["error"].as<double>(),2);
      }
      if (config["as"]) {
        meas[5] = config["as"]["mean"].as<double>();
        var[5]  = pow(config["as"]["error"].as<double>(),2);
      }

      double mu;
      if (config["mu"]) {
        mu = config["mu"].as<double>();
      }
            
      // create FCN function  
      MfitterFcn theFCN(meas, var, mu);

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

