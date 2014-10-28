#include <tt.hpp>
std::complex<long double>
tt::mygl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[29], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=pow(SW,-1);
    aryuttGL[3]=pow(MMH,-1);
    aryuttGL[4]=pow(MMW,-1);
    aryuttGL[5]=pow(MMt,-1);
    aryuttGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuttGL[7]=Tsil::I2(0,0,MMt,mu2);
    aryuttGL[8]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[9]=Tsil::A(MMH,mu2);
    aryuttGL[10]=Tsil::A(MMt,mu2);
    aryuttGL[11]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryuttGL[12]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuttGL[13]=Tsil::Aeps(MMH,mu2);
    aryuttGL[14]=Tsil::Aeps(MMt,mu2);
    aryuttGL[15]=prot0ttHt->M(0);
    aryuttGL[16]=prottH0H->Vxzuv(0);
    aryuttGL[17]=prot0ttHt->Tuxv(0);
    aryuttGL[18]=protWt000->Tyzv(0);
    aryuttGL[19]=protHtt0t->M(0);
    aryuttGL[20]=prot00t00->M(0);
    aryuttGL[21]=prot000t0->M(0);
   aryuttGL[22]=pow(Pi,2);
   aryuttGL[22]= - 17 - 5./2.*aryuttGL[22];
   aryuttGL[23]=1 + 1./2.*aryuttGL[12];
   aryuttGL[23]=aryuttGL[12]*aryuttGL[23];
   aryuttGL[24]=aryuttGL[3]*aryuttGL[13];
   aryuttGL[25]=aryuttGL[10]*aryuttGL[3];
   aryuttGL[26]= - 8*aryuttGL[15] - 8*aryuttGL[19] - aryuttGL[21] - 
   aryuttGL[20];
   aryuttGL[26]=1./3.*aryuttGL[26] + 32*aryuttGL[3];
   aryuttGL[26]=MMt*aryuttGL[26];
   aryuttGL[22]=aryuttGL[26] + 4*aryuttGL[25] + 13./3.*aryuttGL[24] + 
   52./3.*aryuttGL[8] + 5./6.*aryuttGL[23] + 13./3.*aryuttGL[17] - 8./3.
   *aryuttGL[11] + 1./6.*aryuttGL[22] + 4*aryuttGL[18];
   aryuttGL[22]=MMt*aryuttGL[22];
   aryuttGL[23]= - 1 - 7./2.*aryuttGL[8];
   aryuttGL[23]=aryuttGL[10]*aryuttGL[23];
   aryuttGL[23]=aryuttGL[23] - 1./2.*aryuttGL[9] - 5*aryuttGL[14] + 1./
   2.*aryuttGL[6];
   aryuttGL[24]=pow(aryuttGL[10],2);
   aryuttGL[25]=aryuttGL[5]*aryuttGL[24];
   aryuttGL[23]=1./3.*aryuttGL[23] + aryuttGL[25];
   aryuttGL[23]=aryuttGL[5]*aryuttGL[23];
   aryuttGL[25]=4*aryuttGL[11] - 13 - 4*aryuttGL[18];
   aryuttGL[26]=aryuttGL[15] - 8./3.*aryuttGL[16] + aryuttGL[19];
   aryuttGL[26]=MMt*aryuttGL[26];
   aryuttGL[27]=aryuttGL[8] + 1 + aryuttGL[17];
   aryuttGL[27]=aryuttGL[5]*aryuttGL[27];
   aryuttGL[27]=1./2.*aryuttGL[27] - aryuttGL[15] + 4*aryuttGL[16] - 
   aryuttGL[19];
   aryuttGL[27]=MMH*aryuttGL[27];
   aryuttGL[23]=1./3.*aryuttGL[27] + 2*aryuttGL[26] + aryuttGL[23] - 17.
   /3.*aryuttGL[8] + 1./3.*aryuttGL[25] - 2*aryuttGL[17];
   aryuttGL[23]=MMH*aryuttGL[23];
   aryuttGL[25]= - aryuttGL[3]*aryuttGL[9];
   aryuttGL[26]= - aryuttGL[10]*aryuttGL[3];
   aryuttGL[25]=36*aryuttGL[26] + 13./3.*aryuttGL[25] + 14./3.*
   aryuttGL[8] + 23./2. + 4./3.*aryuttGL[12];
   aryuttGL[25]=aryuttGL[10]*aryuttGL[25];
   aryuttGL[26]= - aryuttGL[7] + 5*aryuttGL[14];
   aryuttGL[27]= - aryuttGL[8]*aryuttGL[9];
   aryuttGL[28]= - aryuttGL[3]*pow(aryuttGL[9],2);
   aryuttGL[24]= - aryuttGL[5]*aryuttGL[24];
   aryuttGL[22]=aryuttGL[23] + aryuttGL[22] + 2*aryuttGL[24] + 
   aryuttGL[25] + 4./3.*aryuttGL[28] + 2./3.*aryuttGL[27] + 7./2.*
   aryuttGL[9] - 11./3.*aryuttGL[6] + 13./6.*aryuttGL[26] + 3*
   aryuttGL[13];

      yuttGLret = aryuttGL[22]*aryuttGL[4]*pow(aryuttGL[2],2)*
      aryuttGL[1];
      return yuttGLret;
}
