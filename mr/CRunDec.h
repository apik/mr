/*
  CRunDec.h

  Header file for CRunDec.cpp

  Author: Barbara Schmidt
*/

/*
License:

CRunDec is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef CRUNDEC_H_
#define CRUNDEC_H_

// Default return statement:
#define RETURN return(0);
// The following might be useful for Windows:
//#define RETURN   system("PAUSE"); exit(1);


// Numerical values for input parameters:
#define asMz 0.1183
#define Mz   91.18
#define Mt   173.2
#define Mb   4.8
#define Mc   1.5
#define muc  1.279
#define mub  4.163
#define Mtau 1.777

/* // Some constants: */
#define cf 4./3.
#define ca 3.
#define tr 1./2.
#define B4 -1.762800087073770864061897634679818807215137274389016762629478603776
#define A4 0.5174790616738993863307581618988629456223774751413792582443193479770
#define A5 0.5084005792422687074591088492585899413195411256648216487244977963526
#define Pi M_PI
#define Zeta2 (Pi*Pi)/6.
#define Zeta3 1.20205690315959428539973816151144999076498629234049888179227155534
#define Zeta4 (Pi*Pi*Pi*Pi)/90.
#define Zeta5 1.03692775514336992633136548645703416805708091950191281197419267790


// Struct for triple {nf, Mth, Muth}:
struct TriplenfMmu{
      int nf;
      double Mth;
      double muth;
};
  
struct AsmMS{
      double Asexact;
      double mMSexact;
};

// class declaration of CRunDec:
class CRunDec
{
  private: 
  // Aux. constants for implicit Runhe-Kutta-Procedure:
  static const double a2=0.2, a3=0.3, a4=0.6, a5=1., a6=0.875;
       
  static const double b21=0.2, b31=3./40., b32=9./40., b41=0.3, b42=-0.9, 
                      b43=6./5.;
  static const double b51=-11./54., b52=2.5, b53=-70./27., b54=35./27.;
  static const double b61=1631./55296., b62=175./512., b63=575./13824.;
  static const double b64=44275./110592., b65=253./4096.;
       
  static const double c1=37./378., c2=0., c3=250./621., c4=125./594., c5=0.;
  static const double c6= 512./1771.;
       
  static const double dc1=37./378.-2825./27648., dc2=0.-0., 
                      dc3=250./621.-18575./48384.;
  static const double dc4=125./594.-13525./55296., dc5=0.-277./14336., 
                      dc6=512./1771.-0.25;
  
  // Coefficients for diff. equations:
  double Beta[4], B[4], Gamma[4], C[4], Nf;
  
  // Define constants (if not already done with constructor):
  void SetConstants(int n);
  
  // R.h.s. of diff. equations:
  friend double fSetdydx(CRunDec S, double A,int nl);
  
  friend double fSetdydxa1(CRunDec S,double x, double A);
  friend double fSetdydxM1(CRunDec S,double A, double M);

  friend double fSetdydxa2(CRunDec S,double x, double A);
  friend double fSetdydxM2(CRunDec S,double A, double M);

  friend double fSetdydxa3(CRunDec S,double x, double A);
  friend double fSetdydxM3(CRunDec S,double A, double M);
  
  friend double fSetdydxa4(CRunDec S,double x, double A);
  friend double fSetdydxM4(CRunDec S,double A, double M);
  
  // Additional aux. functions:
  int Abbruch(void);
  double fSetAsL(double Lambda, double Mu, int nl, double AlphaS);
  double fSetcx(double x, int nl);
  double fOsFromMs1(double mu, double M);
  double fOsFromMs2(double mu, double M, double nl);
  double fOsFromMs3(double mu, double M, double nl);

  /* Made public and static to use separately from CRunDec */
  /* A.Pikelner                                            */
 public:
  static double fMsFromOs1(double mu, double M);
  static double fMsFromOs2(double mu, double M, double nl);
  static double fMsFromOs3(double mu, double M, double nl);
 private:
  double fZmM(double n);
  double fZmInvM(double n);
  double fDelta(double mOS,double mq[]);
  double fMsFromRi1(void);
  double fMsFromRi2(void);
  double fMsFromRi3(void);
  double fMumFromOs1(void);
  double fMumFromOs2(void);
  double fMumFromOs3(void);
  double fRiFromMs(double alpha, double nl);
  double fMsFromRi(double alpha, double nl);
  double fHelpmOS2mMSit(double MS,double mOS, double mq[], double asmu,
                        double mu, int nl);
  double fas5to6os(double alpha,double mass, double mu,double nlq, double nl);
  double fas6to5os(double alpha,double mass, double mu,double nlq, double nl);
  double fmq5to6os(double A,double mass,double mu,double nlq,double nl);
  double fmq6to5os(double A,double mass,double mu,double nlq,double nl);
  
  double fRungeKuttaImpl(double &x, double y,double &htry, int nl, 
                      double (*f)(CRunDec,double, int));
  double fRKSchritt(double x,double y,double h,double &yerr,
                    double (*f)(CRunDec, double ,double));
     
  public:
  // constructor:
  CRunDec();
  CRunDec(int);
  
  // Arrays and structs to store data:
  double mq[4];
  TriplenfMmu nfMmu[4];
  AsmMS AM;
  
  // Function to obtain current number of active flavours:
  int GetNf();
  // Function to set number of active flavours:
  void SetNf(int nf);
  
  // Functions for the running of alpha_s:
  double LamExpl(double asmu, double mu, int nloops);
  double LamImpl(double asmu, double mu,int nloops);
  double AlphasLam(double Lambda, double mu, int nloops);
  double AlphasExact(double asmu0, double mu0, double mu1, int nloops);

  // Funktions for the runnung of mq and the
  // various mass definitions:
  double mMS2mMS(double mu0, double asmu0, double asmu1, int nloops);
  double mMS2mOS(double MS, double mq[4], double asmu, double mu, int nloops);
  double mOS2mMS(double mOS, double mq[], double asmu, double mu, int nloops);
  double mMS2mSI(double mMS, double asmu, double mu, int nloops);
  double mRI2mMS(double mRI, double asmu, int nloops);
  double mMS2mRGI(double mMS, double asmu, int nloops);
  double mRGI2mMS(double mRGI, double asmu, int nloops);
  double mOS2mSI(double mOS, double mq[], double asM, int nloops);
  double mOS2mMSrun(double mOS, double mq[], double asmu, double mu, int nloops); 
  double mMS2mOSrun(double mMS, double mq[], double asmu, double mu, int nloops); 
  double mMS2mRI(double mMS, double asmu, int nloops); 
  double mOS2mMSit(double mOS, double mq[], double asmu, double mu, int nloops); 
  double mMS2mRGImod(double mMS, double asmu, int nloops);
  // Solve coupled differential equations for alpha_s and mq:
  AsmMS AsmMSrunexact(double mmu, double asmu0, double mu0, double mu1, 
                      int nloops);

  // Decoupling relations:
  double DecAsDownOS(double asmu, double massth, double muth, int nloops);
  double DecAsUpOS(double asmu, double massth, double muth, int nloops);
  double DecMqUpOS(double mq, double asmu, double massth, double muth, int nloops);
  double DecMqDownOS(double mq, double asmu, double massth, double muth, int nloops);

  // Running and decoupling:
  double AlL2AlH(double asl,double mu1,TriplenfMmu decpar[],double mu2, int nloops);
  double AlH2AlL(double ash,double mu1,TriplenfMmu decpar[],double mu2, int nloops);
  double mL2mH(double mql,double asl,double mu1,TriplenfMmu decpar[],double mu2,
               int nloops);
  double mH2mL(double mqh,double ash,double mu1,TriplenfMmu decpar[],double mu2,
               int nloops);
  
  // Overload functions:
  double LamExpl(double asmu, double mu, int nf, int nloops);
  double LamImpl(double asmu, double mu,int nf,int nloops);
  double AlphasLam(double Lambda, double mu,int nf, int nloops);
  double AlphasExact(double asmu0, double mu0, double mu1, int nf,int nloops);
  double mMS2mMS(double mu0, double asmu1, double asmu0,int nf, int nloops);
  AsmMS AsmMSrunexact(double mmu, double asmu0, double mu0, double mu1,
                       int nf, int nloops);
  double mMS2mOS(double MS, double mq[], double asmu, double mu,int nf, int nloops);
  double mOS2mMS(double mOS, double mq[], double asmu, double mu,int nf,int nloops);
  double mMS2mSI(double mMS, double asmu, double mu,int nf, int nloops);
  double mRI2mMS(double mRI, double asmu,int nf, int nloops);
  double mMS2mRGI(double mMS, double asmu,int nf, int nloops);
  double mRGI2mMS(double mRGI, double asmu,int nf, int nloops);
  double mOS2mSI(double mOS, double mq[], double asM,int nf, int nloops);
  double mOS2mMSrun(double mOS, double mq[], double asmu, double mu,int nf, 
                    int nloops); 
  double mMS2mOSrun(double mMS, double mq[], double asmu, double mu,int nf, 
                    int nloops); 
  double mMS2mRI(double mMS, double asmu,int nf, int nloops); 
  double mOS2mMSit(double mOS, double mq[],double asmu,double mu,int nf,int nloops);
  double mMS2mRGImod(double mMS, double asmu,int nf, int nloops);
  double DecAsDownOS(double asmu, double massth, double muth,int nf, int nloops);
  double DecAsUpOS(double asmu, double massth, double muth, int nf, int nloops);
  double DecMqUpOS(double mq, double asmu, double massth, double muth, int nf,
                   int nloops);
  double DecMqDownOS(double mq, double asmu, double massth, double muth,int nf,
                     int nloops);
  
};

#endif
