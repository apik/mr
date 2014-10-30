//
// MR - 2-loop matching and 3-loop Running, including full 2-loop EW corrections
// Copyright (C) 2014 Andrey Pikelner <pikelner@theor.jinr.ru>
//
// This file is part of MR.
//
// MR is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MR is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MR.  If not, see <http://www.gnu.org/licenses/>.
//

/*
  CRunDec.cpp, v1.1

  Author: Barbara Schmidt

  For documentation see

  CRunDec: a C++ package for running and decoupling of the
  strong coupling and quark masses

  by

  Barbara Schmidt, Matthias Steinhauser

  Comput.Phys.Commun. 183 (2012) 1845-1848
  arXiv:1201.6149 
  SFB/CPP-12-03
  TTP12-02

  See also:
  [RunDec] K.~G.~Chetyrkin, J.~H.~Kuhn and M.~Steinhauser,
  ``RunDec: A Mathematica package for running and decoupling of the strong
  coupling and quark masses,''
  Comput.\ Phys.\ Commun.\  {\bf 133} (2000) 43
  [arXiv:hep-ph/0004189].

  May 2013: minor modifications to avoid "-Wunused-parameter"
            during compilation with g++ -Wall -Wextra
  Sep 2013: bug fix in fRungeKuttaImpl (thanks to Stephen Jones)
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

#include <iostream>
#include <cmath>
#include "CRunDec.h"

using namespace std;

// Default constructor:
CRunDec::CRunDec(){
     for(int i=0; i<4; i++){
       mq[i]=0.;
       nfMmu[i].Mth=0.;
       nfMmu[i].muth=0.;
       nfMmu[i].nf=0;
     }
     AM.Asexact=0.;
     AM.mMSexact=0.;
}

// Constructor called with number for active flavours nf
CRunDec::CRunDec(int n){
     double nf=(double) n;
     Nf=nf;
     Beta[0]= (double)0.25*(11.0 - 2.*nf/3.0);
     Beta[1]= (double)(102.0 - 38.*nf/3.0)/16.0;
     Beta[2]= (double)(0.5*2857. - 5033.*nf/18.0 + 325.*nf*nf/54.0)/64.0;
     Beta[3]= (double)(149753./6.0 + 3564.*Zeta3 + 
              (-1078361./162.0 - 6508.*Zeta3/27.0)*nf + 
              (50065./162.0 + 6472.*Zeta3/81.0)*nf*nf + 
              1093.*nf*nf*nf/729.0)/256.0;

     Gamma[0]=(double)1.0;
     Gamma[1]=(double)(202./3.-20.*nf/9.)/16.;  
     Gamma[2]=(double)(1249. + (-2216./27. - 160.*Zeta3/3.)*nf-
              140.*nf*nf/81.)/64.;
     Gamma[3]=(double)(4603055./162. + 135680.*Zeta3/27. - 8800.*Zeta5 +
              (-91723./27. - 34192.*Zeta3/9. + 
              880.*Zeta4 + 18400.*Zeta5/9.)*nf +
              (5242./243. + 800.*Zeta3/9. - 160.*Zeta4/3.)*nf*nf +
              (-332./243. + 64.*Zeta3/27.)*nf*nf*nf)/256.;
      
     for(int i=0; i<4; i++) {
       B[i]=Beta[i]/Beta[0];
       C[i]=Gamma[i]/Beta[0];
       mq[i]=0.;
       nfMmu[i].Mth=0.;
       nfMmu[i].muth=0.;
       nfMmu[i].nf=0;
     } 
     AM.Asexact=0.;
     AM.mMSexact=0.;      
}

// Define constants (if not already done in constructor)
void CRunDec::SetConstants(int n){
     double nf=(double) n;
     if(Nf!=nf){
     }
     Nf=nf;
     Beta[0]= (double)0.25*(11.0 - 2.*nf/3.0);
     Beta[1]= (double)(102.0 - 38.*nf/3.0)/16.0;
     Beta[2]= (double)(0.5*2857. - 5033.*nf/18.0 + 325.*nf*nf/54.0)/64.0;
     Beta[3]= (double)(149753./6.0 + 3564.*Zeta3 + 
              (-1078361./162.0 - 6508.*Zeta3/27.0)*nf + 
              (50065./162.0 + 6472.*Zeta3/81.0)*nf*nf + 
              1093.*nf*nf*nf/729.0)/256.0;
     
     Gamma[0]=(double)1.0;
     Gamma[1]=(double)(202./3.-20.*nf/9.)/16.; 
     Gamma[2]=(double)(1249. + (-2216./27. - 160.*Zeta3/3.)*nf-
              140.*nf*nf/81.)/64.;
     Gamma[3]=(double)(4603055./162. + 135680.*Zeta3/27. - 8800.*Zeta5 +
              (-91723./27. - 34192.*Zeta3/9. + 
              880.*Zeta4 + 18400.*Zeta5/9.)*nf +
              (5242./243. + 800.*Zeta3/9. - 160.*Zeta4/3.)*nf*nf +
              (-332./243. + 64.*Zeta3/27.)*nf*nf*nf)/256.;
      
     for(int i=0; i<4; i++) {
       B[i]=Beta[i]/Beta[0];
       C[i]=Gamma[i]/Beta[0];
     }                          
}

// Function int CRunDec::GetNf()
// Returns number of active flavours.
int CRunDec::GetNf(){
     return (int)Nf;
}

// Function void CRunDec::SetNf(int nf)
// Set the private component Nf to the number of active flavours.
void CRunDec::SetNf(int nf){
     this->SetConstants(nf);
}

// Aux. function to exit function.
int CRunDec::Abbruch(void){
     RETURN
}

// Function double CRunDec::LamExpl(double AlphaS, double Mu, int nl)
// Compute \Lambda using eq.(4) of [RunDec].
double CRunDec::LamExpl(double AlphaS, double Mu, int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double A=AlphaS/Pi;
     double sum[4];
     sum[0]= 1./(A*Beta[0]);
     sum[1]= (B[1]*log(A))/Beta[0]
              + (B[1]/Beta[0])*log(Beta[0]);
     sum[2]= (B[2]*A - B[1]*B[1]*A )/Beta[0];
     sum[3]= (0.5*B[3]*A*A-B[1]*B[2]*A*A+ 0.5*B[1]*B[1]*B[1]*A*A)/Beta[0];

     double LogM2L2=0.0;
     for(int i=1; i<=nl; i++){
       LogM2L2+=sum[i-1];
     }
               
     double Lambda= Mu*exp(-0.5*LogM2L2);
     return Lambda;
}

// Function double CRunDec::AlphasLam(double Lambda, double Mu, int nl)
// Compute \alpha_s using eq.(5) of [RunDec]. 
double CRunDec::AlphasLam(double Lambda, double Mu, int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(Mu/Lambda<1.5){
       cout<<"WARNING: the ratio \\mu/\\lambda = "<< Mu/Lambda
             <<" is very small!"<<endl; 
       RETURN  
     }
     double L=log(Mu*Mu/(Lambda*Lambda));
     double h=1/(L*Beta[0]);
     double c=log(L);
     double sum[4];
     sum[0]= h;
     sum[1]= -h*h*B[1]*c;
     sum[2]= + h*h*h*(B[1]*B[1]*(c*c-c-1)+B[2]);
     sum[3]= h*h*h*h*(B[1]*B[1]*B[1]*(-c*c*c+2.5*c*c+2*c-0.5) 
                   -3*B[1]*B[2]*c +0.5*B[3]);
     double a=0.0;
     for(int i=1; i<=nl; i++){
       a+=sum[i-1];
     }
     return a*Pi;
}

// Eq.(5) rewritten in a form suitable to determine zero.
double CRunDec::fSetAsL(double Lambda, double Mu, int nl, double AlphaS){
     double L=log(Mu*Mu/(Lambda*Lambda));
     double h=1/(L*Beta[0]);
     double c=log(L);
     double sum[4];
     sum[0]= h;
     sum[1]= -h*h*B[1]*c;
     sum[2]= + h*h*h*(B[1]*B[1]*(c*c-c-1)+B[2]);
     sum[3]= h*h*h*h*(B[1]*B[1]*B[1]*(-c*c*c+2.5*c*c+2*c-0.5) 
                   -3*B[1]*B[2]*c +0.5*B[3]);
     double Add=0.0;
     for(int i=1; i<=nl; i++) Add+=sum[i-1];
     return (Add-(AlphaS/Pi));
}

// Function double CRunDec::LamImpl(double AlphaS, double Mu,int nl)
// Compute \Lambda using eq.(5) of [RunDec]. 
double CRunDec::LamImpl(double AlphaS, double Mu,int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double epsilonX= 1e-8;
     double Lambda0=LamExpl(AlphaS,Mu,nl);
     double x0=Lambda0 - 0.2*Lambda0;
     double x1=Lambda0 + 0.2*Lambda0;
     double f0= this->fSetAsL(x0,Mu,nl,AlphaS);
     double f1= this->fSetAsL(x1,Mu,nl,AlphaS);
     if(f0*f1>0){
       cout<<"WARNING: No root can be calculatet!"<<endl;
       RETURN
     }
     double xTest;
     double fTest;
     do{
       xTest= (x0+x1)/2;
       fTest= fSetAsL(xTest,Mu,nl,AlphaS);
       if(f0*fTest<0){x1= xTest;}
       else {x0= xTest;}
     }
     while(abs(x1-x0)>= epsilonX);
     double Lambda=xTest;
     return Lambda;
}

// Right-hand side of differential equation for \alpha_s times 1/\mu^2,
// see eq.(1) of [1].
double fSetdydx(CRunDec S, double A,int nl){ 
     double f=0.0;
     double sum[4];
     double B=A*A;
     sum[0]=-S.Beta[0]*B;    
     sum[1]=-S.Beta[1]*B*A;
     sum[2]=-S.Beta[2]*B*B;
     sum[3]=-S.Beta[3]*B*B*A;    
     for(int i=1; i<= nl; i++) {
       f+=sum[i-1];
     }
     return (f*2);
}

// Implicit Runge-Kutte step (4th order)
// Call with x=\mu^2, y=\alpha_s(\mu)/\pi, step size h, number of loops nl 
// New y value is returned at x+h
double CRunDec::fRungeKuttaImpl(double &x, double y,double &htry, int nl, 
                              double (*f)(CRunDec,double, int)){
     // Precision
     double eps=1e-10;
     double yerr,ytemp,htemp, hnext;
     double h=htry;
     double xnew; // new variable   
     double k1,k2,k3,k4,k5,k6;
     for(;;){
       k1=h*f(*this,y,nl);
       k2=h*f(*this,y+b21*k1,nl);
       k3=h*f(*this,y+b31*k1+b32*k2,nl);
       k4=h*f(*this,y+b41*k1+b42*k2+b43*k3,nl);
       k5=h*f(*this,y+b51*k1+b52*k2+b53*k3+b54*k4,nl);
       k6=h*f(*this,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5,nl);
       // y value at x+h as a sum of the previous value and the
       // correspondingly weighted function evaluations
       ytemp= y+ c1*k1+ c2*k2+ c3*k3+ c4*k4+ c5*k5+ c6*k6;
       // Estimate of uncertainty
       yerr=dc1*k1 + dc2*k2 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6;
       double err=0.;
       err=fmax(err,fabs(yerr/eps));
       
       // Uncertainty too big? -> Discard result and reduce step size
       if(err>1.){      
         htemp=0.9*h*pow(err,-0.25);
         if(h>=0.){h=fmax(htemp,0.1*h);}
         else{h=fmin(htemp,0.1*h);}
         xnew=x+h;  // modification to previous code
         //decide whether reduced stepsize is still big enough 
         //(in order to prevent a closed loop)
         if(xnew==x){cout<<"stepsize too small"<<endl; RETURN} 
         continue;          
       }
       else{
         if(err>1.89e-4){
         hnext=0.9*h*pow(err,-0.2);
         }
         // Uncertainty OK? -> take y value, increase h
         else{
           hnext=5.*h;
         }
         x+=h;
         
         y=ytemp;
         htry=hnext;
         break;
       }       
     }
     return y; 
}

// Function: double CRunDec::AlphasExact(double AlphaS0, double Mu0, 
//                           double MuEnd, int nl)
// Compute \alpha_s using eq.(1) of [RunDec]
double CRunDec::AlphasExact(double AlphaS0, double Mu0, double MuEnd, int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double Lambda=LamExpl(AlphaS0,Mu0,nl);
     if(MuEnd/Lambda<1.5){
       cout<<"WARNING: the ratio \\mu/\\lambda = "<< MuEnd/Lambda
         <<" is very small!"<<endl;
       RETURN
     }
     double x,y;
     x=log(Mu0);
     y=AlphaS0/Pi;
     double h;
     
     if(Mu0<MuEnd){
       h=1e-4;
       while(x<(log(MuEnd))){
         y=this->fRungeKuttaImpl(x,y,h,nl,fSetdydx);
         if(x+h>=log(MuEnd)){
           h=log(MuEnd)-x;
         }
       }
     return(y*Pi);
     }
     else{h=-1e-4;}
     while(x>(log(MuEnd))){
       y=this->fRungeKuttaImpl(x,y,h,nl,fSetdydx);
       if(x+h<=log(MuEnd)){
         h=log(MuEnd)-x;
       }
     }
     return(y*Pi);
} 

// Eq.(10) of [RunDec]
double CRunDec::fSetcx(double x, int nl){
     double sum[4];
     sum[0]=1;
     sum[1]=(C[1]-B[1]*C[0])*x;
     sum[2]=0.5*((C[1]-B[1]*C[0])*(C[1]-B[1]*C[0]) + C[2] - B[1]*C[1]+
           B[1]*B[1]*C[0] - B[2]*C[0])*x*x;
     sum[3]=((C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])*(C[1]-B[1]*C[0])/6. +
           0.5*(C[1]-B[1]*C[0])*(C[2]-B[1]*C[1]+B[1]*B[1]*C[0]-B[2]*C[0])+
           (C[3]-B[1]*C[2]+B[1]*B[1]*C[1]-B[2]*C[1]-B[1]*B[1]*B[1]*C[0] +
            2.*B[1]*B[2]*C[0] - B[3]*C[0])/3.)*x*x*x;
     double erg=0.0;
     for(int i=1; i<=nl; i++){
       erg+=sum[i-1];
     }
     return (pow(x,C[0])*erg);             
}

// Function double CRunDec::mMS2mMS(double Mu0, double AlphaS0, 
//                          double AlphaSEnd, int nl)
// Compute m_q(\mu) using eqs.(9) and (10) of [RunDec]
double CRunDec::mMS2mMS(double Mu0, double AlphaS0, double AlphaSEnd, int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return Mu0;
     }
     double cAlphaS0= this->fSetcx(AlphaS0/Pi, nl);
     double cAlphaSEnd= this->fSetcx(AlphaSEnd/Pi, nl);
     return Mu0*cAlphaSEnd/cAlphaS0;        
}

// Aux. functions (r.h.s of diff. eqs. for alpha_s and m_q)
double fSetdydxM1(CRunDec S,double A, double M){
     return (M*(S.Gamma[0])/(S.Beta[0]*A));
}

double fSetdydxa1(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A));      
}

double fSetdydxM2(CRunDec S,double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A)/(S.Beta[0]*A+S.Beta[1]*A*A));   
}

double fSetdydxa2(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A));      
}

double fSetdydxM3(CRunDec S,double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A+S.Gamma[2]*A*A)/
           (S.Beta[0]*A+S.Beta[1]*A*A+S.Beta[2]*A*A*A));     
}

double fSetdydxa3(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A+S.Beta[2]*A*A*A*A));      
}

double fSetdydxM4(CRunDec S,double A, double M){
     return (M*(S.Gamma[0]+S.Gamma[1]*A+S.Gamma[2]*A*A+S.Gamma[3]*A*A*A)/
           (S.Beta[0]*A+S.Beta[1]*A*A+S.Beta[2]*A*A*A+S.Beta[3]*A*A*A*A)); 
}

double fSetdydxa4(CRunDec S,double x, double A){
     if (x == 1.) x=1.;
     return (-2.*(S.Beta[0]*A*A+S.Beta[1]*A*A*A+S.Beta[2]*A*A*A*A+
                S.Beta[3]*A*A*A*A*A));      
} 

// Runge-Kutta step for implicit procedure
double CRunDec::fRKSchritt(double x,double y,double h,double &yerr,
                    double (*f)(CRunDec, double ,double)){
     double k1,k2,k3,k4,k5,k6; 
     k1=h*f(*this,x,y);
     k2=h*f(*this,x+a2*h,y+b21*k1);
     k3=h*f(*this,x+a3*h,y+b31*k1+b32*k2);
     k4=h*f(*this,x+a4*h,y+b41*k1+b42*k2+b43*k3);
     k5=h*f(*this,x+a5*h,y+b51*k1+b52*k2+b53*k3+b54*k4);
     k6=h*f(*this,x+a6*h,y+b61*k1+b62*k2+b63*k3+b64*k4+b65*k5);

     yerr=(dc1*k1 + dc2*k2 + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6);
     return (y+ c1*k1+ c2*k2+ c3*k3+ c4*k4+ c5*k5+ c6*k6);      
}  

// Function AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
//                         double MuEnd, int nl)
// Compute \alpha_s and m_q solving diff. eqs. simultaneously.
AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
                             double MuEnd, int nl){
     AsmMS Erg;
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       this->Abbruch();
     }
     if(nl==0){
       Erg.Asexact=AlphaS0;
       Erg.mMSexact=mMu;
       return Erg;
     }
     double yerr0,ytemp0,xnew,h=1e-3;
     float eps=1e-15;
     float errmax;
     double x0=log(Mu0);
     double y0=AlphaS0/Pi;
     double xEnd=log(MuEnd);
     double yscal0=abs(x0)+abs(h*y0);
     
     double (*falpha)(CRunDec, double ,double);
     double (*fmMS)(CRunDec, double ,double);
     
     if(nl==1){
       falpha=fSetdydxa1;
       fmMS=fSetdydxM1;
     }
     if(nl==2){
       falpha=fSetdydxa2;
       fmMS=fSetdydxM2;
     } 
     if(nl==3){
       falpha=fSetdydxa3;
       fmMS=fSetdydxM3;
     }
     if(nl==4){
       falpha=fSetdydxa4;
       fmMS=fSetdydxM4;
     }
     
     if(Mu0<MuEnd){
       h=1e-2;
       while(x0<xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,falpha);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;  
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h>=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.Asexact=y0*Pi;   
       
       x0=AlphaS0/Pi;
       xEnd=y0;
       y0=mMu;
       yscal0=abs(x0)+abs(h*y0);
       eps=1e-10;
       h=1e-3;

       while(x0>xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,fmMS);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h<=xEnd){
           h=xEnd-x0;
         }//if
       }//while  
       Erg.mMSexact=y0;
       return Erg;
     }//if
     else{
       h=-1e-2;
       while(x0>xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,falpha);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h<=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.Asexact=y0*Pi; 

       x0=AlphaS0/Pi;     
       xEnd=y0;
       y0=mMu;
       yscal0=abs(x0)+abs(h*y0);
       eps=1e-10;
       h=1e-3;
       while(x0<xEnd){
         for(;;){
           ytemp0=fRKSchritt(x0,y0,h,yerr0,fmMS);
           errmax=0.;
           errmax=fmax(errmax,fabs((float)yerr0/yscal0));
           errmax/=eps;
           if(errmax>1){
             h*=0.9;
             xnew=x0+h;
             if(xnew==x0){cout<<"stepsize too small!"<<endl;}
             continue;
           } //if
           else{
             x0+=h;
             y0=ytemp0;
             if(errmax>1.89e-4){h=0.9*h*pow((double)errmax,-0.2);}
             else{h=5.*h;}
             break;
           }//else
         }//for
         if(x0+h>=xEnd){
           h=xEnd-x0;
         }//if
       }//while
       Erg.mMSexact=y0;
       return Erg;
     }//else
}   

// Coefficients of eq.(13) of [RunDec]
double CRunDec::fMsFromOs1(double mu, double M){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (-cf - (3.*cf*lmM)/4.);
     return erg;
}

double CRunDec::fMsFromOs2(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) ((-1111.*ca*cf)/384. + (7.*cf*cf)/128. -
     (185.*ca*cf*lmM)/96. + (21.*cf*cf*lmM)/32. - (11.*ca*cf*lmM*lmM)/32. +
     (9.*cf*cf*lmM*lmM)/32. + (143.*cf*tr)/96. + (13.*cf*lmM*tr)/24. + 
     (cf*lmM*lmM*tr)/8. +
     (71.*cf*nl*tr)/96. + (13.*cf*lmM*nl*tr)/24. + (cf*lmM*lmM*nl*tr)/8. +
     (ca*cf*Zeta2)/2 - (15.*cf*cf*Zeta2)/8. - (3.*ca*cf*log(2)*Zeta2)/2. + 
     3.*cf*cf*log(2)*Zeta2 -
     cf*tr*Zeta2 + (cf*nl*tr*Zeta2)/2. + (3.*ca*cf*Zeta3)/8. - 
     (3.*cf*cf*Zeta3)/4.);
     return erg;
}

double CRunDec::fMsFromOs3(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (lmM*lmM*(-2341.*ca*ca*cf + 1962.*ca*cf*cf - 243.*cf*cf*cf 
     + 1492.*ca*cf*tr -
      468.*cf*cf*tr + 1492.*ca*cf*nl*tr - 468.*cf*cf*nl*tr - 208.*cf*tr*tr -
      416.*cf*nl*tr*tr - 208.*cf*nl*nl*tr*tr))/1152. +
     (lmM*lmM*lmM*(-242.*ca*ca*cf + 297.*ca*cf*cf - 81.*cf*cf*cf + 
     176.*ca*cf*tr - 108.*cf*cf*tr + 176.*ca*cf*nl*tr - 108.*cf*cf*nl*tr - 
     32.*cf*tr*tr - 64.*cf*nl*tr*tr - 32.*cf*nl*nl*tr*tr))/1152. +
     (lmM*(-105944.*ca*ca*cf + 52317.*ca*cf*cf - 13203.*cf*cf*cf + 
     74624.*ca*cf*tr -
     5436.*cf*cf*tr + 55616.*ca*cf*nl*tr + 2340.*cf*cf*nl*tr - 
     12608.*cf*tr*tr -
     18304.*cf*nl*tr*tr - 5696.*cf*nl*nl*tr*tr + 12672.*ca*ca*cf*Zeta2 -
     52704.*ca*cf*cf*Zeta2 + 19440.*cf*cf*cf*Zeta2 - 
     38016.*ca*ca*cf*log(2)*Zeta2 +
     91584.*ca*cf*cf*log(2)*Zeta2 - 31104.*cf*cf*cf*log(2)*Zeta2 - 
     29952.*ca*cf*tr*Zeta2 +
     27648.*cf*cf*tr*Zeta2 + 13824.*ca*cf*log(2)*tr*Zeta2 - 
     27648.*cf*cf*log(2)*tr*Zeta2 +
     8064.*ca*cf*nl*tr*Zeta2 + 12096.*cf*cf*nl*tr*Zeta2 + 
     13824.*ca*cf*log(2)*nl*tr*Zeta2 -
     27648.*cf*cf*log(2)*nl*tr*Zeta2 + 9216.*cf*tr*tr*Zeta2 + 
     4608.*cf*nl*tr*tr*Zeta2 -
     4608.*cf*nl*nl*tr*tr*Zeta2 + 9504.*ca*ca*cf*Zeta3 - 
     22896.*ca*cf*cf*Zeta3 +
     7776.*cf*cf*cf*Zeta3 + 6912.*ca*cf*tr*Zeta3 - 3456.*cf*cf*tr*Zeta3 +
     6912.*ca*cf*nl*tr*Zeta3 - 3456.*cf*cf*nl*tr*Zeta3))/13824.;
     return erg;
}

//z[3,m](M) according to eq.(15)
double CRunDec::fZmM(double nl){
     double erg;
     erg= (double) -9478333./93312. + 55.*log(2)*log(2)*log(2)*log(2)/162. +
            (-644201./6480. + 587.*log(2)/27. + 44.*log(2)*log(2)/27.)*Zeta2 -
            61.*Zeta3/27. + 3475*Zeta4/432. + 1439.*Zeta2*Zeta3/72. -
            1975.*Zeta5/216. + 220.*A4/27. + nl*(246643./23328. - 
            log(2)*log(2)*log(2)*log(2)/81. +(967./108. + 22.*log(2)/27. -
            4.*log(2)*log(2)/27.)*Zeta2 + 241.*Zeta3/72. - 305.*Zeta4/108. -
            8.*A4/27.) + nl*nl*(-2353./23328. - 13.*Zeta2/54 - 7.*Zeta3/54.);
     return erg;
}

// Function: double CRunDec::mOS2mMS(double mOS, double mq[], double asmu,
//                           double Mu,int nl)
double CRunDec::mOS2mMS(double mOS, double mq[], double asmu,double Mu,int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     double sum[4];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fMsFromOs1(Mu, mOS))/Pi;
     sum[2]=asmu*asmu*((this-> fMsFromOs2(Mu, mOS, Nf-1))
           -4.*(this->fDelta(mOS,mq)/3.))/(Pi*Pi); 
     sum[3]=asmu*asmu*asmu*(this-> fMsFromOs3(Mu, mOS,Nf-1)+
     this->fZmM(Nf-1))/(Pi*Pi*Pi);
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return mOS*erg;       
}

// Function: double CRunDec::mMS2mSI(double mMS, double asmu, double mu, int nl)
double CRunDec::mMS2mSI(double mMS, double asmu, double mu, int nl){   
  double epsilonX = 1e-8;
  AsmMS asmq;
  asmq.Asexact  = asmu;
  asmq.mMSexact = mMS;
  for (;;) {
    double mbold = asmq.mMSexact;
    asmq = AsmMSrunexact(mMS, asmu, mu, mbold, nl);
    if (abs(asmq.mMSexact - mbold) < epsilonX) break;
  }
  return asmq.mMSexact;
}

// Coefficients of eq.(18) of [RunDec]
double CRunDec::fMsFromRi1(void){
     return (double) -4./3.;
}
double CRunDec::fMsFromRi2(void){
     double erg= (double) -995./72. + 19.*Zeta3/6. + 89.*Nf/144.;
     return erg;       
}
double CRunDec::fMsFromRi3(void){
     double erg= (double) -6663911./41472. + 408007.*Zeta3/6912. -185.*Zeta5/36.
                + (118325./7776. + 5*Zeta4/12. - 617.*Zeta3/216.)*Nf +
                (-4459./23328. - Zeta3/54.)*Nf*Nf;
     return erg;                    
}

// Function: double CRunDec::mRI2mMS(double mRI, double asmu, int nl)
double CRunDec::mRI2mMS(double mRI, double asmu, int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double sum[4];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fMsFromRi1())/Pi;
     sum[2]=asmu*asmu*(this-> fMsFromRi2())/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fMsFromRi3())/(Pi*Pi*Pi);
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
        erg+=sum[i];
       }
     }
     return mRI*erg;         
}

// Function: double CRunDec::mMS2mRGI(double mMS, double asmu, int nl)
double CRunDec::mMS2mRGI(double mMS, double asmu, int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return (double) mMS;
     }
     else{    
       double cAsmu= this->fSetcx(asmu/Pi, nl); 
       return (double) mMS/cAsmu;
     }
}

// Function: double CRunDec::mRGI2mMS(double mRGI, double asmu, int nl)
double CRunDec::mRGI2mMS(double mRGI, double asmu, int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN   
     }  
     if(nl==0){
       return (double) mRGI;
     }
     double cAsmu= this->fSetcx(asmu/Pi, nl); 
     return (double) mRGI*cAsmu;             
}

// Coefficients of eq.(17) of [RunDec]
double CRunDec::fOsFromMs1(double mu, double M){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) (cf + (3.*cf*lmM)/4.);
     return erg;
}

double CRunDec::fOsFromMs2(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;
     erg= (double) ((1111.*ca*cf)/384. - (71.*cf*cf)/128. -
     (143.*cf*tr)/96. - (71.*cf*nl*tr)/96. +
     lmM*((185.*ca*cf)/96. - (9.*cf*cf)/32. - (13.*cf*tr)/24. - 
     (13.*cf*nl*tr)/24.) +
     lmM*lmM*((11.*ca*cf)/32. + (9.*cf*cf)/32. - (cf*tr)/8. - (cf*nl*tr)/8.) -
     (ca*cf*Zeta2)/2. + (15.*cf*cf*Zeta2)/8. + (3.*ca*cf*log(2)*Zeta2)/2. 
     - 3.*cf*cf*log(2)*Zeta2 +
     cf*tr*Zeta2 - (cf*nl*tr*Zeta2)/2. - (3.*ca*cf*Zeta3)/8. + 
     (3.*cf*cf*Zeta3)/4.);
     return erg;
}

double CRunDec::fOsFromMs3(double mu, double M, double nl){
     double lmM=log((mu*mu)/(M*M));
     double erg;

     erg= (double) (lmM*lmM*lmM*((121.*ca*ca*cf)/576. + (33.*ca*cf*cf)/128. +
     (9.*cf*cf*cf)/128. - (11.*ca*cf*tr)/72. - (3.*cf*cf*tr)/32. - 
     (11.*ca*cf*nl*tr)/72. -
     (3.*cf*cf*nl*tr)/32. + (cf*tr*tr)/36. + (cf*nl*tr*tr)/18. +
     (cf*nl*nl*tr*tr)/36.) + lmM*lmM*((2341.*ca*ca*cf)/1152. 
     + (21.*ca*cf*cf)/64. -
     (63.*cf*cf*cf)/128. - (373.*ca*cf*tr)/288. - (3.*cf*cf*tr)/32. -
     (373.*ca*cf*nl*tr)/288. - (3.*cf*cf*nl*tr)/32. + (13.*cf*tr*tr)/72. +
     (13.*cf*nl*tr*tr)/36. + (13.*cf*nl*nl*tr*tr)/72.) +
     lmM*((13243.*ca*ca*cf)/1728. - (4219.*ca*cf*cf)/1536. + 
     (495.*cf*cf*cf)/512. -
     (583.*ca*cf*tr)/108. - (307.*cf*cf*tr)/384. - (869.*ca*cf*nl*tr)/216. -
     (91.*cf*cf*nl*tr)/384. + (197.*cf*tr*tr)/216. + (143.*cf*nl*tr*tr)/108. +
     (89.*cf*nl*nl*tr*tr)/216. - (11.*ca*ca*cf*Zeta2)/12. + 
     (49.*ca*cf*cf*Zeta2)/16. +
     (45.*cf*cf*cf*Zeta2)/32. + (11.*ca*ca*cf*log(2)*Zeta2)/4. - 
     (35.*ca*cf*cf*log(2)*Zeta2)/8. -
     (9.*cf*cf*cf*log(2)*Zeta2)/4. + (13.*ca*cf*tr*Zeta2)/6. - 
     (cf*cf*tr*Zeta2)/2. -
     ca*cf*log(2)*tr*Zeta2 + 2.*cf*cf*log(2)*tr*Zeta2 - 
     (7.*ca*cf*nl*tr*Zeta2)/12. -
     (13.*cf*cf*nl*tr*Zeta2)/8. - ca*cf*log(2)*nl*tr*Zeta2 + 
     2.*cf*cf*log(2)*nl*tr*Zeta2 -
     (2.*cf*tr*tr*Zeta2)/3. - (cf*nl*tr*tr*Zeta2)/3. + 
     (cf*nl*nl*tr*tr*Zeta2)/3. -
     (11.*ca*ca*cf*Zeta3)/16. + (35.*ca*cf*cf*Zeta3)/32. + 
     (9.*cf*cf*cf*Zeta3)/16. -
     (ca*cf*tr*Zeta3)/2. + (cf*cf*tr*Zeta3)/4. - (ca*cf*nl*tr*Zeta3)/2. +
     (cf*cf*nl*tr*Zeta3)/4.));
     return erg;
}

// Compute \Delta using eq.(14) of [RunDec]
double CRunDec::fDelta(double mOS,double mq[]){
     double erg=0.0;
     for (int i=0; i<4; i++){
       double x= (double) mq[i]/mOS;
       if(x>1.||x<0.){
         cout<<"\\Delta "<<x<<" IS CALLED; THE FUNCTION IS NOT IMPELEMENTED "<<
               "FOR ARGUMENTS OUTSIDE THE INTERVAL [0,1]."<<endl;
         RETURN
       }
       erg+= Pi*Pi*x/8. - 0.597*x*x + 0.230*x*x*x;
     }
     return erg;
}       
// z_m^inv
double CRunDec::fZmInvM(double nl){
     double erg;
     erg=(8481925./93312. + 
       (137.*nl)/216. + (652841.*Pi*Pi)/38880. - (nl*Pi*Pi)/27. - 
       (695.*Pi*Pi*Pi*Pi)/7776. - (575.*Pi*Pi*log(2))/162. - 
       (22.*Pi*Pi*log(2)*log(2))/81. - 
       (55.*log(2)*log(2)*log(2)*log(2))/162. - (220.*A4)/27. - 
       nl*nl*(-2353./23328. - (13.*Pi*Pi)/324. - (7.*Zeta3)/54.) + 
       (58.*Zeta3)/27. - 
       (1439.*Pi*Pi*Zeta3)/432. - nl*(246643./23328. + (967.*Pi*Pi)/648. - 
       (61.*Pi*Pi*Pi*Pi)/1944. + (11.*Pi*Pi*log(2))/81. - 
       (2.*Pi*Pi*log(2)*log(2))/81. - 
       log(2)*log(2)*log(2)*log(2)/81. - (8.*A4)/27. + 
       (241.*Zeta3)/72.) + 
       (1975.*Zeta5)/216.);
     return erg;
}

// Function: double CRunDec::mMS2mOS(double MS, double mq[], double asmu,
//                           double mu,int nl)
double CRunDec::mMS2mOS(double MS, double mq[], double asmu,double mu,int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double sum[4];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fOsFromMs1(mu, MS))/Pi;
     sum[2]=asmu*asmu*((this-> fOsFromMs2(mu, MS, Nf-1))
           +4.*(this->fDelta(MS,mq))/3.)/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fOsFromMs3(mu, MS,Nf-1)+
           this->fZmInvM(Nf-1))/(Pi*Pi*Pi);     
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{ 
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return MS*erg;
}    

// Coefficients of eq.(16) of [RunDec]
double CRunDec::fMumFromOs1(void){
     return (double) -cf;
}
  
double CRunDec::fMumFromOs2(void){
     double erg;
     erg= (double) ((-1111.*ca*cf)/384. + (199.*cf*cf)/128. + (143.*cf*tr)/96. +
     (71.*cf*(Nf-1)*tr)/96. + (ca*cf*Zeta2)/2. - (15.*cf*cf*Zeta2)/8. - 
     (3.*ca*cf*log(2)*Zeta2)/2. +
     3.*cf*cf*log(2)*Zeta2 - cf*tr*Zeta2 + (cf*(Nf-1)*tr*Zeta2)/2. + 
     (3.*ca*cf*Zeta3)/8. -
     (3.*cf*cf*Zeta3)/4.);     
     return erg;      
}

// z_m^SI
double CRunDec::fMumFromOs3(void){
     double erg;
     erg= (double) -7172965./93312. - 
	 (293.*(Nf-1))/216. - (618281.*Pi*Pi)/38880. - ((Nf-1)*Pi*Pi)/9. + 
     (695.*Pi*Pi*Pi*Pi)/7776. + (623.*Pi*Pi*log(2))/162. + 
     (22.*Pi*Pi*log(2)*log(2))/81. + 
     (55.*log(2)*log(2)*log(2)*log(2))/162. + (220.*A4)/27. + 
     (Nf-1)*(Nf-1)*(-2353./23328. - (13.*Pi*Pi)/324. - (7.*Zeta3)/54.) - 
     (70.*Zeta3)/27. + 
     (1439.*Pi*Pi*Zeta3)/432. + (Nf-1)*(246643./23328. + (967.*Pi*Pi)/648. - 
     (61.*Pi*Pi*Pi*Pi)/1944. + (11.*Pi*Pi*log(2))/81. - 
     (2*Pi*Pi*log(2)*log(2))/81. - 
     log(2)*log(2)*log(2)*log(2)/81. - (8.*A4)/27. + (241.*Zeta3)/72.) - 
     (1975.*Zeta5)/216.;
     return erg;   
}

// Function: double CRunDec::mOS2mSI(double mOS, double mq[], double asM, 
//                           int nl)
double CRunDec::mOS2mSI(double mOS, double mq[], double asM, int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double sum[4];
     sum[0]=(double) 1.;
     sum[1]=asM*(this ->fMumFromOs1())/Pi;
     sum[2]=asM*asM*(this-> fMumFromOs2()-4.*fDelta(mOS,mq)/3.)/(Pi*Pi);
     sum[3]=asM*asM*asM*(this-> fMumFromOs3())/(Pi*Pi*Pi);       
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return mOS*erg;       
       
}

// Function: double CRunDec::mOS2mMSrun(double mOS, double mq[], double asmu, 
//                           double mu, int nl)
double CRunDec::mOS2mMSrun(double mOS, double mq[], double asmu, double mu,
			   int nl){
     double asM=0.0;
     asM= this-> AlphasExact(asmu, mu, mOS, nl);
     double mum= this-> mOS2mSI(mOS, mq, asM, nl);
     double asmum= this-> AlphasExact(asmu, mu, mum, nl);
     double newM= this->mMS2mMS(mum, asmum, asmu, nl);
     return newM;       
}

// Function: double CRunDec::mMS2mOSrun(double mMS, double mq[], double asmu, 
//                           double mu, int nl)
double CRunDec::mMS2mOSrun(double mMS, double mq[], double asmu, double mu,
			   int nl){
     double mNeu = mMS2mSI(mMS, asmu, mu, nl);
     double asmNeu = AlphasExact(asmu, mu, mNeu, nl);
     return mMS2mOS(mNeu, mq, asmNeu, mNeu, nl);   
}

// Coefficients of eq.(19) of [RunDec]
double CRunDec::fRiFromMs(double alpha, double nl){
     double sum[4];
     sum[0]= 1.;
     sum[1]= (4.*alpha)/3.;
     sum[2]= alpha*alpha*((1123./72. - (89.*Nf)/144. - (19.*Zeta3)/6.));
     sum[3]= alpha*alpha*alpha*(6663911./41472. - (118325.*Nf)/7776. +
            (4459.*Nf*Nf)/23328. + 
            (4.*(1123./72. - (89.*Nf)/144. - (19.*Zeta3)/6.))/3. - 
            (408007.*Zeta3)/6912. + 
            (617.*Nf*Zeta3)/216. + (Nf*Nf*Zeta3)/54. - (4.*(-995./72. + 
            (89.*Nf)/144. + (19.*Zeta3)/6.))/3. - 
            (5.*Nf*Zeta4)/12. + (185.*Zeta5)/36.);
     double erg=0.0;
     if(nl==0){
       erg=1.;
     }
     else{
       erg=0.;
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return (erg);        
}

// Function: double CRunDec::mMS2mRI(double mMS, double asmu, int nl)
double CRunDec::mMS2mRI(double mMS, double asmu, int nl){
     if(nl<0||nl>3){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     return (double) mMS*(this->fRiFromMs((asmu/Pi), nl)); 
}

// Coefficients needed for the transformation of mOS to mMSit
double CRunDec::fHelpmOS2mMSit(double MS, double mOS, double mq[], double asmu,
                                double mu,int nl){
     double sum[4];
     sum[0]=(double) 1.;
     sum[1]=asmu*(this ->fOsFromMs1(mu, MS))/Pi;
     sum[2]=asmu*asmu*((this-> fOsFromMs2(mu, MS, Nf-1))
           +4.*(this->fDelta(mOS,mq))/3.)/(Pi*Pi);
     sum[3]=asmu*asmu*asmu*(this-> fOsFromMs3(mu, MS,Nf-1)+
           +this->fZmInvM(Nf-1))/(Pi*Pi*Pi);
     double erg=0.0;
     if(nl==0){
       erg=1;
     }
     else{ 
       for(int i=0; i<=nl; i++){
         erg+=sum[i];
       }
     }
     return erg;
}

// Function: double CRunDec::mOS2mMSit(double mOS, double mq[], double asmu, 
//                          double mu, int nl)
double CRunDec::mOS2mMSit(double mOS, double mq[], double asmu, double mu,
                           int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN
     }
     double epsilonX= 1e-8;
     double x0=mOS-0.1*mOS; 
     double x1=mOS+0.1*mOS;
     double f0= (x0*(this->fHelpmOS2mMSit(x0,mOS, mq, asmu, mu, nl))-mOS);
     double f1= (x1*(this->fHelpmOS2mMSit(x1,mOS, mq, asmu, mu, nl))-mOS);
     if(f0*f1>0){
       cout<<"WARNING: No root can be calculatet!"<<endl;
       RETURN
     }
     double xTest;
     double fTest;
     do{
       xTest= (x0+x1)/2.;
       fTest= (xTest*(this->fHelpmOS2mMSit(xTest,mOS, mq, asmu, mu, nl))-mOS);
       if(f0*fTest<=0){x1= xTest;}
       else {x0= xTest;}
     }
     while(abs(x1-x0)>= epsilonX);
       double mNeu=xTest;
     return mNeu;               
}

// Function: double CRunDec::mMS2mRGImod(double mMS, double asmu, int nl)
// See 'mMS2mRGI' but 'Alphas/Pi -> AlphaS*2*Beta0/Pi'
double CRunDec::mMS2mRGImod(double mMS, double asmu, int nl){
     if(nl<0||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     if(nl==0){
       return (double) mMS;
     }
     double bet0= 11./4. - Nf/6.;     
     double cAsmu= this->fSetcx((2*bet0*asmu)/Pi, nl); 
     return (double) mMS/cAsmu;     
}

// Coefficients of eq.(22) of [RunDec]
double CRunDec::fas5to6os(double A,double mass,double mu,double nlq,double nl){
     double log2 = log(2.);
     double lmM  = log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=A*(-lmM)/6.;
     sum[2]=A*A*(-7./24. - (19.*lmM)/24. + (lmM*lmM)/36.);
     sum[3]=A*A*A*(-58933./124416. - (8521.*lmM)/1728. - (131.*lmM*lmM)/576. -
                     (lmM*lmM*lmM)/216. +
       nlq*(2479./31104. + (409.*lmM)/1728. + Zeta2/9.) - (2.*Zeta2)/3. - 
       (2.*Zeta2*log(2))/9. -
       (80507.*Zeta3)/27648.);
     sum[4]=
       + pow(A,4)*(
		   -39.499192691910871413 + (2965705*A4)/54432. - (121*A5)/36. - 
       (7693*pow(lmM,2))/1152. - (8371*pow(lmM,3))/10368. + 
       pow(lmM,4)/1296. - (697121*Zeta2)/19440. - (49*log2*Zeta2)/54. - 
       (243892631*Zeta4)/8.70912e6 + (605*Zeta2*Zeta4)/2688. + 
       (330575*Zeta5)/41472. + (587*Zeta2*log2)/81. - 
       (2057*Zeta4*log2)/576. - (2699593*Zeta2*pow(log2,2))/217728. - 
       (121*Zeta2*pow(log2,3))/432. + 
       lmM*(-26.38581988383059 - (29*Zeta2)/9. - (29*log2*Zeta2)/27. - 
          (2439119*Zeta3)/165888.) + 
       pow(nlq,2)*(-0.09432401513203018 - (493*pow(lmM,2))/20736. + 
          lmM*(-0.008996699245541839 - Zeta2/27.) - (13*Zeta2)/162. - 
          (19*Zeta3)/1728.) + (1439*Zeta2*Zeta3)/216. + 
       nlq*(2.375194240826475 + (173*A4)/5184. + (6661*pow(lmM,2))/10368. + 
          (107*pow(lmM,3))/1728. + (557*Zeta2)/162. - 
          (697709*Zeta4)/165888. + (115*Zeta5)/576. + (22*Zeta2*log2)/81. - 
	    (1709*Zeta2*pow(log2,2))/20736. + 
          (173*pow(log2,4))/124416. + (4756441*Zeta3)/995328. + 
          lmM*(2.975080911351166 + (41*Zeta2)/54. + (2*log2*Zeta2)/27. + 
             (132283*Zeta3)/82944.))
		   );

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsDownOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsDownOS(double als, double massth, double muth, int nl){
       if(nl<1||nl>5){
         cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
         RETURN 
       }
       double erg=(this->fas5to6os(als/Pi, massth, muth, Nf,nl));
       return als*erg;
}

// Coefficients of eq.(25) of [RunDec]
double CRunDec::fas6to5os(double A,double mass,double mu,double nlq,double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[5];
     sum[0]=1.;
     sum[1]=A*(lmM)/6.;
     sum[2]=A*A*(7./24. + (19.*lmM)/24. + (lmM*lmM)/36.);
     sum[3]=A*A*A*(58933./124416. + (2.*Zeta2)/3. + (2.*Zeta2*log(2))/9. + 
                    (80507.*Zeta3)/27648. + (8941.*lmM)/1728. +
                    (511.*lmM*lmM)/576. + (lmM*lmM*lmM)/216. +
                    nlq*(-2479./31104. - Zeta2/9. - (409.*lmM)/1728));
     sum[4]=A*A*A*A*
       (39.754401025244204746 + (47039*pow(lmM,2))/3456 + 
	(14149*pow(lmM,3))/10368 + pow(lmM,4)/1296 - 
	(2965705*A4)/54432 + (121  *  A5  )/36 + 
	(697121*Zeta2)/19440 + (49*log(2)*Zeta2)/54 + (243892631*Zeta4)/8709120 - 
	(605*Zeta2*Zeta4)/2688 - (330575*Zeta5)/41472 - (587*Zeta2*log(2))/81 + 
	(2057*Zeta4*log(2))/576 + (2699593*Zeta2*pow(log(2),2))/217728 + 
	(121*Zeta2*pow(log(2),3))/432 + 
	nlq*(-1773073/746496 - (9115*pow(lmM,2))/10368 - 
	     (107*pow(lmM,3))/1728 - (173*A4)/5184 
	     - (557*Zeta2)/162 + (697709*Zeta4)/165888 - 
	     (115*Zeta5)/576 - (22*Zeta2*log(2))/81 + 
	     (1709*Zeta2*pow(log(2),2))/20736 - 
	     (173*pow(log(2),4))/124416 + 
	     lmM*(-1140191/373248 - (47*Zeta2)/54 - 
		  (2*log(2)*Zeta2)/27 - (132283*Zeta3)/82944) 
	     - (4756441*Zeta3)/995328) + 
	nlq*nlq*(140825/1492992 + (493*pow(lmM,2))/20736 + 
		 lmM*(1679/186624 + Zeta2/27) + 
		 (13*Zeta2)/162 + (19*Zeta3)/1728) - (1439*Zeta2*Zeta3)/216 + 
	lmM*(21084715/746496 + (35*Zeta2)/9 + (35*log(2)*Zeta2)/27 + 
	     (2922161*Zeta3)/165888));

     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecAsUpOS(double als, double massth, double muth, int nl){
     if(nl<1||nl>5){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl <<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fas6to5os(als/Pi, massth, muth, Nf,nl));
     return als*erg;
}

// Coefficients of eq.(33) of [RunDec]
double CRunDec::fmq6to5os(double A,double mass,double mu,double nlq,double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[4];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(-89./432. + (5.*lmM)/36. - (lmM*lmM)/12.);
     sum[3]=A*A*A*(-1871./2916. + 407.*Zeta3/864. - 5.*Zeta4/4. + B4/36. -
                     299.*lmM/2592. + 5.*Zeta3*lmM/6. - 299*lmM*lmM/432. - 
                     35.*lmM*lmM*lmM/216. + nlq*(-1327./11664. + 2.*Zeta3/27. +
                     53.*lmM/432. + lmM*lmM*lmM/108.));
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqUpOS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN 
     }
     double erg=(this->fmq6to5os(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Coefficients of eq.(30) of [RunDec]
double CRunDec::fmq5to6os(double A,double mass,double mu,double nlq,double nl){
     double lmM=log((mu*mu)/(mass*mass));
     if(nl==1){
       return 1;
     }
     double sum[4];
     sum[0]=1.;
     sum[1]=0.;
     sum[2]=A*A*(89./432. - (5.*lmM)/36. + (lmM*lmM)/12.);
     sum[3]=A*A*A*(+1871./2916. - 407.*Zeta3/864. + 5.*Zeta4/4. - B4/36. +
                     121.*lmM/2592. - 5.*Zeta3*lmM/6. + 319*lmM*lmM/432. + 
                     29.*lmM*lmM*lmM/216. + nlq*(1327./11664. - 2.*Zeta3/27. -
                     53.*lmM/432. - lmM*lmM*lmM/108.));
     double erg=0.0;
     for(int i=1;i<=nl;i++){
       erg+=sum[i-1];
     }
     return erg;
}

// Function double CRunDec::DecAsUpOS(double als, double massth, double muth, 
//                          int nl)
double CRunDec::DecMqDownOS(double mq, double als, double massth, double muth, 
                           int nl){
     if(nl<1||nl>4){
       cout<<"PROCEDURE IS NOT IMPLEMENTED FOR "<< nl<<" LOOPS"<<endl;
       RETURN  
     }
     double erg=(this->fmq5to6os(als/Pi, massth, muth, Nf,nl));
     return mq*erg;
}

// Function double CRunDec::AlL2AlH(double asl,double mu1,TriplenfMmu decpar[],
//                          double mu2, int nl)
double CRunDec::AlL2AlH(double asl,double mu1,TriplenfMmu decpar[],double mu2, 
                         int nl){
     int n=0;
     int help;
     double help2;
     double asini=asl;
     double muini=mu1;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf>decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=3-n+2;i<=3;i++){
       if(decpar[i].nf-decpar[i-1].nf!=1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT."<<endl;
         RETURN  
       }
     }
     double erg1,erg2;
     int i;
     for(i=3-n+1;i<=3;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf-1,nl);
       erg2= DecAsUpOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       muini=decpar[i].muth; 
     }
     double alpha= AlphasExact(asini,muini,mu2,decpar[i-1].nf,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return alpha;
}

// Function double CRunDec::AlH2AlL(double ash,double mu1,TriplenfMmu decpar[],
//                          double mu2, int nl)
double CRunDec::AlH2AlL(double ash,double mu1,TriplenfMmu decpar[],double mu2, 
                         int nl){
     int n=0;
     int help;
     double help2;
     double asini=ash;
     double muini=mu1;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf<decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=1;i<=n-1;i++){
       if(decpar[i].nf-decpar[i-1].nf!=-1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     double erg1,erg2;
     int i;
     for(i=0;i<=n-1;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf,nl);
       erg2= DecAsDownOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       muini=decpar[i].muth; 
     }
     double alpha= AlphasExact(asini,muini,mu2,decpar[i-1].nf-1,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return alpha;
}

// Function double CRunDec::mL2mH(double mql,double asl,double mu1,
//                         TriplenfMmu decpar[], double mu2, int nl)
double CRunDec::mL2mH(double mql,double asl,double mu1,TriplenfMmu decpar[],
                       double mu2, int nl){
     int n=0;
     int help;
     double help2;
     double asini=asl;
     double muini=mu1;
     double mqini=mql;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf>decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=3-n+2;i<=3;i++){
       if(decpar[i].nf-decpar[i-1].nf!=1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     
     double erg1,erg2,erg3, erg4, erg5, erg6;
     int i;
     for(i=3-n+1;i<=3;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf-1,nl);
       erg3= mMS2mMS(mqini,asini,erg1,decpar[i].nf-1,nl);
       erg2= DecAsUpOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       erg4=DecMqUpOS(erg3,erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       asini=erg2;
       mqini=erg4;
       muini=decpar[i].muth; 
     }
     erg5= AlphasExact(asini,muini,mu2,decpar[i-1].nf,nl);
     erg6= mMS2mMS(mqini,asini,erg5,decpar[i-1].nf,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return erg6;
}

// Function double CRunDec::mH2mL(double mqh,double ash,double mu1,
//                         TriplenfMmu decpar[], double mu2, int nl)
double CRunDec::mH2mL(double mqh,double ash,double mu1,TriplenfMmu decpar[],
                      double mu2, int nl){
     int n=0;
     int help;
     double help2;
     double asini=ash;
     double muini=mu1;
     double mqini=mqh;
     for(int i=0; i<4; i++){
       if(decpar[i].nf!=0){
         n+=1;
       }
     }
     int zaehler=4;
     while (zaehler!=0){
       for(int i=0; i<zaehler-1; i++){
         if(decpar[i].nf<decpar[i+1].nf){
           help=decpar[i].nf;
           decpar[i].nf=decpar[i+1].nf;
           decpar[i+1].nf=help;
           help2=decpar[i].Mth;
           decpar[i].Mth=decpar[i+1].Mth;
           decpar[i+1].Mth=help2;
           help2=decpar[i].muth;
           decpar[i].muth=decpar[i+1].muth;
           decpar[i+1].muth=help2;
         }
       }
       zaehler-=1;
     }
     for(int i=1;i<=n-1;i++){
       if(decpar[i].nf-decpar[i-1].nf!=-1){
         cout<<"WARNING: THERE IS A GAP IN NUMBER OF FLAVOURS. EXIT.";
         RETURN  
       }
     }
     double erg1,erg2,erg3,erg4,erg5,erg6;
     int i;
     for( i=0;i<=n-1;i++){
       erg1= AlphasExact(asini, muini, decpar[i].muth,decpar[i].nf,nl);
       erg3= mMS2mMS(mqini,asini,erg1,decpar[i].nf,nl);
       erg2= DecAsDownOS(erg1,decpar[i].Mth,decpar[i].muth,decpar[i].nf-1,nl);
       erg4=DecMqDownOS(erg3,erg1,decpar[i].Mth,decpar[i].muth,
            decpar[i].nf-1,nl);
       asini=erg2;
       mqini=erg4;
       muini=decpar[i].muth; 
     }
     erg5= AlphasExact(asini,muini,mu2,decpar[i-1].nf-1,nl);
     erg6 =mMS2mMS(mqini,asini,erg5,decpar[i-1].nf-1,nl);
     for(int j=0;j<=3;j++){
       decpar[j].nf=0;
       decpar[j].Mth=0.;
       decpar[j].muth=0.;
     }
     return erg6;
}


// In the following the functions above are overloaded w.r.t. to an
// additional argument, n_f (number of active flavours).
// Use SetConstants(nf)
double CRunDec::LamExpl(double AlphaS, double Mu, int nf, int nl){
     SetConstants(nf);
     return (this->LamExpl(AlphaS, Mu, nl));      
}
  
double CRunDec::LamImpl(double AlphaS, double Mu,int nf,int nl){
     SetConstants(nf);
     return (this->LamImpl(AlphaS, Mu, nl));
}
  
double CRunDec::AlphasLam(double Lambda, double Mu,int nf, int nl){
     SetConstants(nf);
     return (this->AlphasLam(Lambda, Mu, nl));
}
  
double CRunDec::AlphasExact(double AlphaS0, double Mu0, double MuEnd, 
                                int nf,int nl){
     SetConstants(nf);
     return (this->AlphasExact(AlphaS0,Mu0,MuEnd,nl));
}
  
double CRunDec::mMS2mMS(double Mu0, double AlphaSEnd, double AlphaS0,int nf, 
                           int nl){
     SetConstants(nf);
     return (this->mMS2mMS(Mu0, AlphaSEnd, AlphaS0, nl));
}
  
AsmMS CRunDec::AsmMSrunexact(double mMu, double AlphaS0, double Mu0, 
                              double MuEnd,int nf, int nl){
     SetConstants(nf);
     return (this->AsmMSrunexact(mMu, AlphaS0, Mu0, MuEnd, nl));
}
  
double CRunDec::mMS2mOS(double MS, double mq[], double asmu,double mu,int nf,
                           int nl){
     SetConstants(nf);
     return (this->mMS2mOS(MS, mq, asmu, mu, nl));
} 
  
double CRunDec::mOS2mMS(double mOS, double mq[],double asmu,double Mu,int nf,
                           int nl){
     SetConstants(nf);
     return (this->mOS2mMS(mOS, mq, asmu, Mu, nl));
}
  
double CRunDec::mMS2mSI(double mMS, double asmu, double mu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mSI(mMS, asmu, mu, nl));
}
  
double CRunDec::mRI2mMS(double mRI, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mRI2mMS(mRI, asmu, nl));
}
  
double CRunDec::mMS2mRGI(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRGI(mMS, asmu, nl));
}
  
double CRunDec::mRGI2mMS(double mRGI, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mRGI2mMS(mRGI, asmu, nl));
}
  
double CRunDec::mOS2mSI(double mOS, double mq[], double asM,int nf, int nl){
     SetConstants(nf);
     return (this->mOS2mSI(mOS, mq, asM, nl));
}
  
double CRunDec::mOS2mMSrun(double mOS, double mq[], double asmu, double mu,
                              int nf, int nl){
     SetConstants(nf);
     return (this->mOS2mMSrun(mOS, mq, asmu, mu, nl));
}
   
double CRunDec::mMS2mOSrun(double mMS, double mq[], double asmu, double mu,
                              int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mOSrun(mMS, mq, asmu, mu, nl));
}
   
double CRunDec::mMS2mRI(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRI(mMS, asmu, nl));
}
  
double CRunDec::mOS2mMSit(double mOS, double mq[], double asmu, double mu,
                             int nf,int nl){
     SetConstants(nf);
     return (this->mOS2mMSit(mOS, mq, asmu, mu, nl));
}
  
double CRunDec::mMS2mRGImod(double mMS, double asmu,int nf, int nl){
     SetConstants(nf);
     return (this->mMS2mRGImod(mMS, asmu, nl));
}       

double CRunDec::DecAsDownOS(double als, double massth, double muth,int nf, 
                               int nl){
     SetConstants(nf);
     return (this->DecAsDownOS(als, massth, muth, nl));
}
  
double CRunDec::DecAsUpOS(double als, double massth, double muth, int nf,
                             int nl){
     SetConstants(nf);
     return (this->DecAsUpOS(als, massth, muth, nl)); 
}
  
double CRunDec::DecMqUpOS(double mq, double als, double massth, double muth,
                             int nf, int nl){
     SetConstants(nf);
     return (this->DecMqUpOS(mq, als, massth, muth, nl));
}
  
double CRunDec::DecMqDownOS(double mq, double als, double massth, 
                               double muth,int nf, int nl){
     SetConstants(nf);
     return (this->DecMqDownOS(mq, als, massth, muth, nl));
}

