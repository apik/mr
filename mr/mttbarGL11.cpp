#include <tt.hpp>
std::complex<long double>
tt::mgl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[29], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=pow(SW,-1);
    armttbarGL[3]=pow(MMH,-1);
    armttbarGL[4]=pow(MMW,-1);
    armttbarGL[5]=pow(MMt,-1);
    armttbarGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    armttbarGL[7]=Tsil::I2(0,0,MMt,mu2);
    armttbarGL[8]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbarGL[9]=Tsil::A(MMH,mu2);
    armttbarGL[10]=Tsil::A(MMt,mu2);
    armttbarGL[11]=Tsil::Beps(MMH,MMt,MMt,mu2);
    armttbarGL[12]=std::real(Tsil::B(0,0,MMt,mu2));
    armttbarGL[13]=Tsil::Aeps(MMH,mu2);
    armttbarGL[14]=Tsil::Aeps(MMt,mu2);
    armttbarGL[15]=prot0ttHt->M(0);
    armttbarGL[16]=prottH0H->Vxzuv(0);
    armttbarGL[17]=prot0ttHt->Tuxv(0);
    armttbarGL[18]=protWt000->Tyzv(0);
    armttbarGL[19]=protHtt0t->M(0);
    armttbarGL[20]=prot00t00->M(0);
    armttbarGL[21]=prot000t0->M(0);
   armttbarGL[22]=pow(Pi,2);
   armttbarGL[22]= - 17 - 5./2.*armttbarGL[22];
   armttbarGL[23]=1 + 1./2.*armttbarGL[12];
   armttbarGL[23]=armttbarGL[12]*armttbarGL[23];
   armttbarGL[24]=armttbarGL[3]*armttbarGL[13];
   armttbarGL[25]=armttbarGL[10]*armttbarGL[3];
   armttbarGL[26]= - 8*armttbarGL[15] - 8*armttbarGL[19] - 
   armttbarGL[21] - armttbarGL[20];
   armttbarGL[26]=1./3.*armttbarGL[26] + 32*armttbarGL[3];
   armttbarGL[26]=MMt*armttbarGL[26];
   armttbarGL[22]=armttbarGL[26] + 4*armttbarGL[25] + 13./3.*
   armttbarGL[24] + 52./3.*armttbarGL[8] + 5./6.*armttbarGL[23] + 13./3.
   *armttbarGL[17] - 8./3.*armttbarGL[11] + 1./6.*armttbarGL[22] + 4*
   armttbarGL[18];
   armttbarGL[22]=MMt*armttbarGL[22];
   armttbarGL[23]= - 1 - 7./2.*armttbarGL[8];
   armttbarGL[23]=armttbarGL[10]*armttbarGL[23];
   armttbarGL[23]=armttbarGL[23] - 1./2.*armttbarGL[9] - 5*
   armttbarGL[14] + 1./2.*armttbarGL[6];
   armttbarGL[24]=pow(armttbarGL[10],2);
   armttbarGL[25]=armttbarGL[5]*armttbarGL[24];
   armttbarGL[23]=1./3.*armttbarGL[23] + armttbarGL[25];
   armttbarGL[23]=armttbarGL[5]*armttbarGL[23];
   armttbarGL[25]=4*armttbarGL[11] - 13 - 4*armttbarGL[18];
   armttbarGL[26]=armttbarGL[15] - 8./3.*armttbarGL[16] + 
   armttbarGL[19];
   armttbarGL[26]=MMt*armttbarGL[26];
   armttbarGL[27]=armttbarGL[8] + 1 + armttbarGL[17];
   armttbarGL[27]=armttbarGL[5]*armttbarGL[27];
   armttbarGL[27]=1./2.*armttbarGL[27] - armttbarGL[15] + 4*
   armttbarGL[16] - armttbarGL[19];
   armttbarGL[27]=MMH*armttbarGL[27];
   armttbarGL[23]=1./3.*armttbarGL[27] + 2*armttbarGL[26] + 
   armttbarGL[23] - 17./3.*armttbarGL[8] + 1./3.*armttbarGL[25] - 2*
   armttbarGL[17];
   armttbarGL[23]=MMH*armttbarGL[23];
   armttbarGL[25]= - armttbarGL[3]*armttbarGL[9];
   armttbarGL[26]= - armttbarGL[10]*armttbarGL[3];
   armttbarGL[25]=36*armttbarGL[26] + 13./3.*armttbarGL[25] + 14./3.*
   armttbarGL[8] + 23./2. + 4./3.*armttbarGL[12];
   armttbarGL[25]=armttbarGL[10]*armttbarGL[25];
   armttbarGL[26]= - armttbarGL[7] + 5*armttbarGL[14];
   armttbarGL[27]= - armttbarGL[8]*armttbarGL[9];
   armttbarGL[28]= - armttbarGL[3]*pow(armttbarGL[9],2);
   armttbarGL[24]= - armttbarGL[5]*armttbarGL[24];
   armttbarGL[22]=armttbarGL[23] + armttbarGL[22] + 2*armttbarGL[24] + 
   armttbarGL[25] + 4./3.*armttbarGL[28] + 2./3.*armttbarGL[27] + 7./2.
   *armttbarGL[9] - 11./3.*armttbarGL[6] + 13./6.*armttbarGL[26] + 3*
   armttbarGL[13];

      mttbarGLret = armttbarGL[22]*armttbarGL[4]*pow(armttbarGL[2],2)*
      armttbarGL[1];
      return mttbarGLret;
}
