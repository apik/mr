#include <tt.hpp>
long double tt::xgl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[30], mttbarGLret;

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
   armttbarGL[22]=armttbarGL[5]*MMH;
   armttbarGL[23]=1./2.*armttbarGL[22];
   armttbarGL[24]=armttbarGL[6] + MMH;
   armttbarGL[24]=armttbarGL[24]*armttbarGL[23];
   armttbarGL[25]=pow(MMH,2);
   armttbarGL[26]=1./2.*armttbarGL[5];
   armttbarGL[26]=armttbarGL[25]*armttbarGL[26];
   armttbarGL[23]=2 - armttbarGL[23];
   armttbarGL[23]=armttbarGL[10]*armttbarGL[23];
   armttbarGL[23]=7*armttbarGL[23] - 2*armttbarGL[9] - 17*MMH + 
   armttbarGL[26];
   armttbarGL[23]=armttbarGL[8]*armttbarGL[23];
   armttbarGL[26]=armttbarGL[15] + armttbarGL[19];
   armttbarGL[27]=4*armttbarGL[16] - armttbarGL[26];
   armttbarGL[27]=MMH*armttbarGL[27];
   armttbarGL[27]=armttbarGL[27] - 13 + 4*armttbarGL[11];
   armttbarGL[27]=MMH*armttbarGL[27];
   armttbarGL[28]=13./2. - armttbarGL[22];
   armttbarGL[28]=armttbarGL[14]*armttbarGL[28];
   armttbarGL[23]=armttbarGL[23] + 5*armttbarGL[28] + armttbarGL[24] - 
   11*armttbarGL[6] + armttbarGL[27];
   armttbarGL[24]=1./3.*armttbarGL[22];
   armttbarGL[22]= - 2 + armttbarGL[22];
   armttbarGL[22]=armttbarGL[10]*armttbarGL[5]*armttbarGL[22];
   armttbarGL[22]=armttbarGL[22] + 4./3.*armttbarGL[12] + 23./2. - 
   armttbarGL[24];
   armttbarGL[22]=armttbarGL[10]*armttbarGL[22];
   armttbarGL[27]=2*MMH;
   armttbarGL[25]=armttbarGL[5]*armttbarGL[25];
   armttbarGL[25]= - armttbarGL[27] + 1./6.*armttbarGL[25];
   armttbarGL[25]=armttbarGL[17]*armttbarGL[25];
   armttbarGL[28]=pow(armttbarGL[9],2);
   armttbarGL[29]= - 13./3.*armttbarGL[9] - 36*armttbarGL[10];
   armttbarGL[29]=armttbarGL[10]*armttbarGL[29];
   armttbarGL[28]= - 4./3.*armttbarGL[28] + armttbarGL[29];
   armttbarGL[28]=armttbarGL[3]*armttbarGL[28];
   armttbarGL[24]=7 - armttbarGL[24];
   armttbarGL[24]=armttbarGL[9]*armttbarGL[24];
   armttbarGL[29]=armttbarGL[18]*MMH;
   armttbarGL[22]= - 13./6.*armttbarGL[7] + armttbarGL[28] + 
   armttbarGL[25] + armttbarGL[22] - 4./3.*armttbarGL[29] + 1./2.*
   armttbarGL[24] + 3*armttbarGL[13] + 1./3.*armttbarGL[23];
   armttbarGL[23]=pow(armttbarGL[2],2);
   armttbarGL[22]=armttbarGL[23]*armttbarGL[22];
   armttbarGL[24]= - 8./3.*armttbarGL[16] + armttbarGL[26];
   armttbarGL[24]=armttbarGL[24]*armttbarGL[27];
   armttbarGL[25]= - 17./2. - 8*armttbarGL[11];
   armttbarGL[27]=1 + 1./2.*armttbarGL[12];
   armttbarGL[27]=armttbarGL[12]*armttbarGL[27];
   armttbarGL[28]=pow(Pi,2);
   armttbarGL[29]=13./3.*armttbarGL[13] + 4*armttbarGL[10];
   armttbarGL[29]=armttbarGL[3]*armttbarGL[29];
   armttbarGL[26]= - armttbarGL[21] - 8*armttbarGL[26] - armttbarGL[20]
   ;
   armttbarGL[26]=32*armttbarGL[3] + 1./3.*armttbarGL[26];
   armttbarGL[26]=MMt*armttbarGL[26];
   armttbarGL[24]=armttbarGL[26] + 52./3.*armttbarGL[8] + 
   armttbarGL[29] + 13./3.*armttbarGL[17] + 4*armttbarGL[18] - 5./12.*
   armttbarGL[28] + 5./6.*armttbarGL[27] + 1./3.*armttbarGL[25] + 
   armttbarGL[24];
   armttbarGL[23]=MMt*armttbarGL[23]*armttbarGL[24];
   armttbarGL[22]=armttbarGL[23] + armttbarGL[22];

      mttbarGLret = armttbarGL[22]*armttbarGL[4]*armttbarGL[1];
      return mttbarGLret.real();
}
