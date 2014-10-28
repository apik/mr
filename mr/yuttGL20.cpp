#include <tt.hpp>
std::complex<long double>
tt::mygl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuttGL[47], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=pow(SW,-1);
    aryuttGL[3]=pow(MMH,-1);
    aryuttGL[4]=pow(MMW,-1);
    aryuttGL[5]=pow(MMt,-1);
    aryuttGL[6]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuttGL[7]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuttGL[8]=Tsil::I2(0,0,MMH,mu2);
    aryuttGL[9]=Tsil::I2(0,0,MMt,mu2);
    aryuttGL[10]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuttGL[11]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[12]=Tsil::A(MMH,mu2);
    aryuttGL[13]=Tsil::A(MMt,mu2);
    aryuttGL[14]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuttGL[15]=std::real(Tsil::B(0,0,MMH,mu2));
    aryuttGL[16]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuttGL[17]=Tsil::B(0,MMH,MMt,mu2);
    aryuttGL[18]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryuttGL[19]=Tsil::Aeps(MMH,mu2);
    aryuttGL[20]=Tsil::Aeps(MMt,mu2);
    aryuttGL[21]=prot0ttHt->Tuxv(0);
    aryuttGL[22]=protHHttH->M(0);
    aryuttGL[23]=protHttHt->M(0);
    aryuttGL[24]=protHHttH->Uxzuv(0);
    aryuttGL[25]=protHHttH->Suxv(0);
    aryuttGL[26]=protHHttH->Uuyxv(0);
    aryuttGL[27]=protWt000->Tyzv(0);
    aryuttGL[28]=protH0tt0->M(0);
    aryuttGL[29]=protH0t00->M(0);
    aryuttGL[30]=prot0Htt0->M(0);
    aryuttGL[31]=prot0H0t0->M(0);
    aryuttGL[32]=prot00ttH->M(0);
    aryuttGL[33]=protH0t00->Uxzuv(0);
    aryuttGL[34]=protH0t00->Uzxyv(0);
    aryuttGL[35]=protH0tt0->Uuyxv(0);
    aryuttGL[36]=protH0t00->Txuv(0);
   aryuttGL[37]=3*aryuttGL[34] - 39./2. + aryuttGL[35];
   aryuttGL[38]= - aryuttGL[8] + aryuttGL[20];
   aryuttGL[38]=aryuttGL[5]*aryuttGL[38];
   aryuttGL[39]=aryuttGL[12]*aryuttGL[5];
   aryuttGL[40]= - 3*aryuttGL[10] - 5 - aryuttGL[15];
   aryuttGL[40]=3*aryuttGL[40] + aryuttGL[11];
   aryuttGL[40]=aryuttGL[11]*aryuttGL[40];
   aryuttGL[41]=aryuttGL[13]*aryuttGL[5];
   aryuttGL[42]= - aryuttGL[17] - 1 - aryuttGL[36];
   aryuttGL[42]=MMH*aryuttGL[5]*aryuttGL[42];
   aryuttGL[37]=1./2.*aryuttGL[42] + 1./2.*aryuttGL[41] + 1./2.*
   aryuttGL[40] + 1./2.*aryuttGL[39] + 1./2.*aryuttGL[38] + 27./4.*
   aryuttGL[10] + 9./4.*aryuttGL[15] + 5./4.*aryuttGL[17] - 21./4.*
   aryuttGL[18] - 23./4.*aryuttGL[21] + 3./4.*aryuttGL[24] + 1./2.*
   aryuttGL[37] + 3*aryuttGL[26];
   aryuttGL[37]=MMH*aryuttGL[37];
   aryuttGL[38]= - aryuttGL[12]*aryuttGL[5];
   aryuttGL[39]= - aryuttGL[13]*aryuttGL[5];
   aryuttGL[40]=3*aryuttGL[15];
   aryuttGL[41]=9*aryuttGL[10];
   aryuttGL[38]=3*aryuttGL[39] - 23*aryuttGL[11] + aryuttGL[38] + 
   aryuttGL[41] + 31 + aryuttGL[40];
   aryuttGL[38]=aryuttGL[13]*aryuttGL[38];
   aryuttGL[39]= - 3*aryuttGL[6] - 5*aryuttGL[8];
   aryuttGL[39]= - 13*aryuttGL[7] + 17./2.*aryuttGL[19] + 1./2.*
   aryuttGL[39] + 31*aryuttGL[20];
   aryuttGL[42]=3./4.*aryuttGL[10] + 3 + 1./4.*aryuttGL[15];
   aryuttGL[42]=aryuttGL[12]*aryuttGL[42];
   aryuttGL[43]= - aryuttGL[11]*aryuttGL[12];
   aryuttGL[37]=aryuttGL[37] + 1./2.*aryuttGL[38] + 1./2.*aryuttGL[43]
    + 1./2.*aryuttGL[39] + 3*aryuttGL[42];
   aryuttGL[37]=MMH*aryuttGL[37];
   aryuttGL[38]=pow(Pi,2);
   aryuttGL[38]= - 473./4. + 3*aryuttGL[38];
   aryuttGL[38]=75./2.*aryuttGL[27] + 1./4.*aryuttGL[38] - 5*
   aryuttGL[33];
   aryuttGL[39]=1 + aryuttGL[16];
   aryuttGL[39]=aryuttGL[16]*aryuttGL[39];
   aryuttGL[42]= - 57 + 11./2.*aryuttGL[16];
   aryuttGL[42]=1./2.*aryuttGL[42] + 15*aryuttGL[14];
   aryuttGL[42]=1./2.*aryuttGL[42] + 3*aryuttGL[11];
   aryuttGL[42]=aryuttGL[11]*aryuttGL[42];
   aryuttGL[38]=aryuttGL[42] - 21./2.*aryuttGL[14] + 3./32.*
   aryuttGL[39] - 3./32.*aryuttGL[17] + 13./8.*aryuttGL[18] - 115./8.*
   aryuttGL[21] + 1./4.*aryuttGL[36] + 1./8.*aryuttGL[38] - 
   aryuttGL[24];
   aryuttGL[39]= - 33 + aryuttGL[16];
   aryuttGL[42]= - 3*aryuttGL[14];
   aryuttGL[39]=1./4.*aryuttGL[39] + aryuttGL[42];
   aryuttGL[39]=aryuttGL[12]*aryuttGL[39];
   aryuttGL[43]=1./2.*aryuttGL[9] - 2*aryuttGL[20];
   aryuttGL[44]=aryuttGL[11]*aryuttGL[12];
   aryuttGL[39]=1./2.*aryuttGL[44] + 1./4.*aryuttGL[39] + 9./4.*
   aryuttGL[7] + 3*aryuttGL[43] - 185./32.*aryuttGL[19];
   aryuttGL[39]=aryuttGL[3]*aryuttGL[39];
   aryuttGL[43]= - 3*aryuttGL[11] - 5./8.*aryuttGL[16] + aryuttGL[14];
   aryuttGL[43]=aryuttGL[3]*aryuttGL[43];
   aryuttGL[44]=pow(aryuttGL[3],2);
   aryuttGL[45]= - aryuttGL[13]*aryuttGL[44];
   aryuttGL[43]=aryuttGL[43] + 3./2.*aryuttGL[45];
   aryuttGL[43]=aryuttGL[13]*aryuttGL[43];
   aryuttGL[45]=5 - aryuttGL[16];
   aryuttGL[46]= - 3 - aryuttGL[14];
   aryuttGL[46]=aryuttGL[11]*aryuttGL[46];
   aryuttGL[45]=aryuttGL[46] + 1./2.*aryuttGL[45] + aryuttGL[14];
   aryuttGL[45]=aryuttGL[3]*aryuttGL[45];
   aryuttGL[45]=aryuttGL[23] + 3*aryuttGL[45];
   aryuttGL[44]= - aryuttGL[13]*aryuttGL[44]*aryuttGL[14];
   aryuttGL[44]=1./2.*aryuttGL[45] + 18*aryuttGL[44];
   aryuttGL[44]=MMt*aryuttGL[44];
   aryuttGL[45]=aryuttGL[31] + aryuttGL[29];
   aryuttGL[45]= - aryuttGL[23] + 1./4.*aryuttGL[45] + 3*aryuttGL[22];
   aryuttGL[45]=MMH*aryuttGL[45];
   aryuttGL[38]=aryuttGL[44] + 1./4.*aryuttGL[45] + 3*aryuttGL[43] + 1./
   4.*aryuttGL[38] + aryuttGL[39];
   aryuttGL[38]=MMt*aryuttGL[38];
   aryuttGL[39]=aryuttGL[41] + 211./2. + aryuttGL[40];
   aryuttGL[39]=1./4.*aryuttGL[39] + aryuttGL[42];
   aryuttGL[39]=1./2.*aryuttGL[39] - aryuttGL[11];
   aryuttGL[39]=aryuttGL[11]*aryuttGL[39];
   aryuttGL[42]=3./2.*aryuttGL[36] - 9*aryuttGL[26] + 5*aryuttGL[27] + 
   5./2.*aryuttGL[33] + 89./2. - 11*aryuttGL[34];
   aryuttGL[42]=9*aryuttGL[14] - 9./2.*aryuttGL[10] - 3./2.*
   aryuttGL[15] - 3./8.*aryuttGL[17] + 35./4.*aryuttGL[18] + 1./2.*
   aryuttGL[42] + 35*aryuttGL[21];
   aryuttGL[43]= - aryuttGL[23] - 9*aryuttGL[22] + aryuttGL[28] + 
   aryuttGL[32] + aryuttGL[30];
   aryuttGL[43]=MMH*aryuttGL[43];
   aryuttGL[39]=1./8.*aryuttGL[43] + 1./4.*aryuttGL[42] + aryuttGL[39];
   aryuttGL[39]=MMH*aryuttGL[39];
   aryuttGL[42]=85./2.*aryuttGL[19] + 17./2.*aryuttGL[25] + 5*
   aryuttGL[9];
   aryuttGL[43]= - 69./2. + aryuttGL[16];
   aryuttGL[43]=1./4.*aryuttGL[43] + aryuttGL[14];
   aryuttGL[43]=aryuttGL[12]*aryuttGL[43];
   aryuttGL[42]=3./2.*aryuttGL[43] + 1./4.*aryuttGL[42] + 3*aryuttGL[7]
   ;
   aryuttGL[40]=11*aryuttGL[11] + aryuttGL[14] + 1./2.*aryuttGL[16] + 
   aryuttGL[41] - 145./16. + aryuttGL[40];
   aryuttGL[41]=aryuttGL[3]*aryuttGL[12];
   aryuttGL[40]=3*aryuttGL[40] + 203./4.*aryuttGL[41];
   aryuttGL[40]=aryuttGL[13]*aryuttGL[40];
   aryuttGL[41]=pow(aryuttGL[12],2);
   aryuttGL[43]=aryuttGL[3]*aryuttGL[41];
   aryuttGL[39]=aryuttGL[39] + 1./2.*aryuttGL[40] + 1./2.*aryuttGL[42]
    + aryuttGL[43];
   aryuttGL[38]=1./4.*aryuttGL[39] + aryuttGL[38];
   aryuttGL[38]=MMt*aryuttGL[38];
   aryuttGL[39]= - 51./32.*aryuttGL[12] + aryuttGL[13];
   aryuttGL[39]=aryuttGL[13]*aryuttGL[39];
   aryuttGL[37]=aryuttGL[38] + 1./16.*aryuttGL[37] + 1./128.*
   aryuttGL[41] + aryuttGL[39];

      yuttGLret = aryuttGL[37]*pow(aryuttGL[4],2)*pow(aryuttGL[2],4)*
      aryuttGL[1];
      return yuttGLret;
}
