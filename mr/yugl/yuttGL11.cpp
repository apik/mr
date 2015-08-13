#include <tt.hpp>
long double tt<OS>::ygl11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[32], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=pow(SW,-1);
    aryuttGL[3]=pow(MMW,-1);
    aryuttGL[4]=pow(MMt,-1);
    aryuttGL[5]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuttGL[6]=Tsil::I2(0,0,MMt,mu2);
    aryuttGL[7]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[8]=Tsil::A(MMH,mu2);
    aryuttGL[9]=Tsil::A(MMt,mu2);
    aryuttGL[10]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryuttGL[11]=pow(MMH,-1);
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
   aryuttGL[22]=aryuttGL[4]*aryuttGL[9];
   aryuttGL[23]=MMH*aryuttGL[4];
   aryuttGL[24]=1./2.*aryuttGL[23] - 17 - 7./2.*aryuttGL[22];
   aryuttGL[24]=MMH*aryuttGL[24];
   aryuttGL[25]=26*MMt + 7*aryuttGL[9] - aryuttGL[8];
   aryuttGL[24]=2*aryuttGL[25] + aryuttGL[24];
   aryuttGL[24]=aryuttGL[7]*aryuttGL[24];
   aryuttGL[25]= - MMt*aryuttGL[20];
   aryuttGL[25]=aryuttGL[25] + 59./4. + 13*aryuttGL[17];
   aryuttGL[25]=MMt*aryuttGL[25];
   aryuttGL[26]=5./4.*aryuttGL[12] + 5./2.;
   aryuttGL[26]=MMt*aryuttGL[26];
   aryuttGL[26]=4*aryuttGL[9] + aryuttGL[26];
   aryuttGL[26]=aryuttGL[12]*aryuttGL[26];
   aryuttGL[27]=pow(MMt,2);
   aryuttGL[28]=aryuttGL[21]*aryuttGL[27];
   aryuttGL[24]= - aryuttGL[25] - aryuttGL[26] + aryuttGL[28] - 
   aryuttGL[24];
   aryuttGL[25]=aryuttGL[5] - aryuttGL[8];
   aryuttGL[26]=1./6.*aryuttGL[4];
   aryuttGL[25]=aryuttGL[26]*aryuttGL[25];
   aryuttGL[26]=1./2.*aryuttGL[4];
   aryuttGL[28]=1 + aryuttGL[17];
   aryuttGL[26]=aryuttGL[28]*aryuttGL[26];
   aryuttGL[26]=4*aryuttGL[16] + aryuttGL[26];
   aryuttGL[28]=1./3.*MMH;
   aryuttGL[26]=aryuttGL[26]*aryuttGL[28];
   aryuttGL[29]=aryuttGL[4]*pow(aryuttGL[9],2);
   aryuttGL[30]= - 1./12.*aryuttGL[9] + aryuttGL[29];
   aryuttGL[30]=aryuttGL[4]*aryuttGL[30];
   aryuttGL[31]=MMt*aryuttGL[16];
   aryuttGL[25]=aryuttGL[26] - 16./3.*aryuttGL[31] + aryuttGL[30] - 53./
   12. - 2*aryuttGL[17] + aryuttGL[25];
   aryuttGL[25]=MMH*aryuttGL[25];
   aryuttGL[26]=2*MMt;
   aryuttGL[30]=aryuttGL[26] - aryuttGL[28];
   aryuttGL[30]=aryuttGL[30]*MMH;
   aryuttGL[27]=aryuttGL[30] - 8./3.*aryuttGL[27];
   aryuttGL[30]=aryuttGL[15] + aryuttGL[19];
   aryuttGL[27]=aryuttGL[27]*aryuttGL[30];
   aryuttGL[30]= - 13*aryuttGL[9] - 4*aryuttGL[8];
   aryuttGL[30]=aryuttGL[11]*aryuttGL[30];
   aryuttGL[22]=1./3.*aryuttGL[30] + 4 - 3./2.*aryuttGL[22];
   aryuttGL[22]=aryuttGL[8]*aryuttGL[22];
   aryuttGL[29]=aryuttGL[29] - aryuttGL[9];
   aryuttGL[23]=13./2. - aryuttGL[23];
   aryuttGL[23]=aryuttGL[14]*aryuttGL[23];
   aryuttGL[30]=MMt*aryuttGL[11];
   aryuttGL[30]=3 + 13./3.*aryuttGL[30];
   aryuttGL[30]=aryuttGL[13]*aryuttGL[30];
   aryuttGL[31]=MMt*pow(Pi,2);
   aryuttGL[26]= - aryuttGL[26] + MMH;
   aryuttGL[26]=aryuttGL[10]*aryuttGL[26];
   aryuttGL[28]=MMt - aryuttGL[28];
   aryuttGL[28]=aryuttGL[18]*aryuttGL[28];
   aryuttGL[22]=4*aryuttGL[28] + 4./3.*aryuttGL[26] - 3./4.*
   aryuttGL[31] + aryuttGL[30] + 5./3.*aryuttGL[23] + aryuttGL[27] + 
   aryuttGL[25] - 11./3.*aryuttGL[5] - 13./6.*aryuttGL[6] + 
   aryuttGL[22] - 8*aryuttGL[29] - 1./3.*aryuttGL[24];

      yuttGLret = aryuttGL[22]*aryuttGL[3]*pow(aryuttGL[2],2)*
      aryuttGL[1];
      return yuttGLret.real();
}
