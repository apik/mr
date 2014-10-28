#include <tt.hpp>
std::complex<long double>
tt::mgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[47], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=pow(SW,-1);
    armttbarGL[3]=pow(MMH,-1);
    armttbarGL[4]=pow(MMW,-1);
    armttbarGL[5]=pow(MMt,-1);
    armttbarGL[6]=Tsil::I2(MMH,MMH,MMH,mu2);
    armttbarGL[7]=Tsil::I2(MMH,MMt,MMt,mu2);
    armttbarGL[8]=Tsil::I2(0,0,MMH,mu2);
    armttbarGL[9]=Tsil::I2(0,0,MMt,mu2);
    armttbarGL[10]=Tsil::B(MMH,MMH,MMH,mu2);
    armttbarGL[11]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbarGL[12]=Tsil::A(MMH,mu2);
    armttbarGL[13]=Tsil::A(MMt,mu2);
    armttbarGL[14]=Tsil::B(MMt,MMt,MMH,mu2);
    armttbarGL[15]=std::real(Tsil::B(0,0,MMH,mu2));
    armttbarGL[16]=std::real(Tsil::B(0,0,MMt,mu2));
    armttbarGL[17]=Tsil::B(0,MMH,MMt,mu2);
    armttbarGL[18]=Tsil::Beps(MMH,MMt,MMt,mu2);
    armttbarGL[19]=Tsil::Aeps(MMH,mu2);
    armttbarGL[20]=Tsil::Aeps(MMt,mu2);
    armttbarGL[21]=prot0ttHt->Tuxv(0);
    armttbarGL[22]=protHHttH->M(0);
    armttbarGL[23]=protHttHt->M(0);
    armttbarGL[24]=protHHttH->Uxzuv(0);
    armttbarGL[25]=protHHttH->Suxv(0);
    armttbarGL[26]=protHHttH->Uuyxv(0);
    armttbarGL[27]=protWt000->Tyzv(0);
    armttbarGL[28]=protH0tt0->M(0);
    armttbarGL[29]=protH0t00->M(0);
    armttbarGL[30]=prot0Htt0->M(0);
    armttbarGL[31]=prot0H0t0->M(0);
    armttbarGL[32]=prot00ttH->M(0);
    armttbarGL[33]=protH0t00->Uxzuv(0);
    armttbarGL[34]=protH0t00->Uzxyv(0);
    armttbarGL[35]=protH0tt0->Uuyxv(0);
    armttbarGL[36]=protH0t00->Txuv(0);
   armttbarGL[37]=3*armttbarGL[34] - 39./2. + armttbarGL[35];
   armttbarGL[38]= - armttbarGL[8] + armttbarGL[20];
   armttbarGL[38]=armttbarGL[5]*armttbarGL[38];
   armttbarGL[39]=armttbarGL[12]*armttbarGL[5];
   armttbarGL[40]= - 3*armttbarGL[10] - 5 - armttbarGL[15];
   armttbarGL[40]=3*armttbarGL[40] + armttbarGL[11];
   armttbarGL[40]=armttbarGL[11]*armttbarGL[40];
   armttbarGL[41]=armttbarGL[13]*armttbarGL[5];
   armttbarGL[42]= - armttbarGL[17] - 1 - armttbarGL[36];
   armttbarGL[42]=MMH*armttbarGL[5]*armttbarGL[42];
   armttbarGL[37]=1./2.*armttbarGL[42] + 1./2.*armttbarGL[41] + 1./2.*
   armttbarGL[40] + 1./2.*armttbarGL[39] + 1./2.*armttbarGL[38] + 27./4.
   *armttbarGL[10] + 9./4.*armttbarGL[15] + 5./4.*armttbarGL[17] - 21./
   4.*armttbarGL[18] - 23./4.*armttbarGL[21] + 3./4.*armttbarGL[24] + 1.
   /2.*armttbarGL[37] + 3*armttbarGL[26];
   armttbarGL[37]=MMH*armttbarGL[37];
   armttbarGL[38]= - armttbarGL[12]*armttbarGL[5];
   armttbarGL[39]= - armttbarGL[13]*armttbarGL[5];
   armttbarGL[40]=3*armttbarGL[15];
   armttbarGL[41]=9*armttbarGL[10];
   armttbarGL[38]=3*armttbarGL[39] - 23*armttbarGL[11] + armttbarGL[38]
    + armttbarGL[41] + 31 + armttbarGL[40];
   armttbarGL[38]=armttbarGL[13]*armttbarGL[38];
   armttbarGL[39]= - 3*armttbarGL[6] - 5*armttbarGL[8];
   armttbarGL[39]= - 13*armttbarGL[7] + 17./2.*armttbarGL[19] + 1./2.*
   armttbarGL[39] + 31*armttbarGL[20];
   armttbarGL[42]=3./4.*armttbarGL[10] + 3 + 1./4.*armttbarGL[15];
   armttbarGL[42]=armttbarGL[12]*armttbarGL[42];
   armttbarGL[43]= - armttbarGL[11]*armttbarGL[12];
   armttbarGL[37]=armttbarGL[37] + 1./2.*armttbarGL[38] + 1./2.*
   armttbarGL[43] + 1./2.*armttbarGL[39] + 3*armttbarGL[42];
   armttbarGL[37]=MMH*armttbarGL[37];
   armttbarGL[38]=pow(Pi,2);
   armttbarGL[38]= - 473./4. + 3*armttbarGL[38];
   armttbarGL[38]=75./2.*armttbarGL[27] + 1./4.*armttbarGL[38] - 5*
   armttbarGL[33];
   armttbarGL[39]=1 + armttbarGL[16];
   armttbarGL[39]=armttbarGL[16]*armttbarGL[39];
   armttbarGL[42]= - 57 + 11./2.*armttbarGL[16];
   armttbarGL[42]=1./2.*armttbarGL[42] + 15*armttbarGL[14];
   armttbarGL[42]=1./2.*armttbarGL[42] + 3*armttbarGL[11];
   armttbarGL[42]=armttbarGL[11]*armttbarGL[42];
   armttbarGL[38]=armttbarGL[42] - 21./2.*armttbarGL[14] + 3./32.*
   armttbarGL[39] - 3./32.*armttbarGL[17] + 13./8.*armttbarGL[18] - 115.
   /8.*armttbarGL[21] + 1./4.*armttbarGL[36] + 1./8.*armttbarGL[38] - 
   armttbarGL[24];
   armttbarGL[39]= - 33 + armttbarGL[16];
   armttbarGL[42]= - 3*armttbarGL[14];
   armttbarGL[39]=1./4.*armttbarGL[39] + armttbarGL[42];
   armttbarGL[39]=armttbarGL[12]*armttbarGL[39];
   armttbarGL[43]=1./2.*armttbarGL[9] - 2*armttbarGL[20];
   armttbarGL[44]=armttbarGL[11]*armttbarGL[12];
   armttbarGL[39]=1./2.*armttbarGL[44] + 1./4.*armttbarGL[39] + 9./4.*
   armttbarGL[7] + 3*armttbarGL[43] - 185./32.*armttbarGL[19];
   armttbarGL[39]=armttbarGL[3]*armttbarGL[39];
   armttbarGL[43]= - 3*armttbarGL[11] - 5./8.*armttbarGL[16] + 
   armttbarGL[14];
   armttbarGL[43]=armttbarGL[3]*armttbarGL[43];
   armttbarGL[44]=pow(armttbarGL[3],2);
   armttbarGL[45]= - armttbarGL[13]*armttbarGL[44];
   armttbarGL[43]=armttbarGL[43] + 3./2.*armttbarGL[45];
   armttbarGL[43]=armttbarGL[13]*armttbarGL[43];
   armttbarGL[45]=5 - armttbarGL[16];
   armttbarGL[46]= - 3 - armttbarGL[14];
   armttbarGL[46]=armttbarGL[11]*armttbarGL[46];
   armttbarGL[45]=armttbarGL[46] + 1./2.*armttbarGL[45] + 
   armttbarGL[14];
   armttbarGL[45]=armttbarGL[3]*armttbarGL[45];
   armttbarGL[45]=armttbarGL[23] + 3*armttbarGL[45];
   armttbarGL[44]= - armttbarGL[13]*armttbarGL[44]*armttbarGL[14];
   armttbarGL[44]=1./2.*armttbarGL[45] + 18*armttbarGL[44];
   armttbarGL[44]=MMt*armttbarGL[44];
   armttbarGL[45]=armttbarGL[31] + armttbarGL[29];
   armttbarGL[45]= - armttbarGL[23] + 1./4.*armttbarGL[45] + 3*
   armttbarGL[22];
   armttbarGL[45]=MMH*armttbarGL[45];
   armttbarGL[38]=armttbarGL[44] + 1./4.*armttbarGL[45] + 3*
   armttbarGL[43] + 1./4.*armttbarGL[38] + armttbarGL[39];
   armttbarGL[38]=MMt*armttbarGL[38];
   armttbarGL[39]=armttbarGL[41] + 211./2. + armttbarGL[40];
   armttbarGL[39]=1./4.*armttbarGL[39] + armttbarGL[42];
   armttbarGL[39]=1./2.*armttbarGL[39] - armttbarGL[11];
   armttbarGL[39]=armttbarGL[11]*armttbarGL[39];
   armttbarGL[42]=3./2.*armttbarGL[36] - 9*armttbarGL[26] + 5*
   armttbarGL[27] + 5./2.*armttbarGL[33] + 89./2. - 11*armttbarGL[34];
   armttbarGL[42]=9*armttbarGL[14] - 9./2.*armttbarGL[10] - 3./2.*
   armttbarGL[15] - 3./8.*armttbarGL[17] + 35./4.*armttbarGL[18] + 1./2.
   *armttbarGL[42] + 35*armttbarGL[21];
   armttbarGL[43]= - armttbarGL[23] - 9*armttbarGL[22] + armttbarGL[28]
    + armttbarGL[32] + armttbarGL[30];
   armttbarGL[43]=MMH*armttbarGL[43];
   armttbarGL[39]=1./8.*armttbarGL[43] + 1./4.*armttbarGL[42] + 
   armttbarGL[39];
   armttbarGL[39]=MMH*armttbarGL[39];
   armttbarGL[42]=85./2.*armttbarGL[19] + 17./2.*armttbarGL[25] + 5*
   armttbarGL[9];
   armttbarGL[43]= - 69./2. + armttbarGL[16];
   armttbarGL[43]=1./4.*armttbarGL[43] + armttbarGL[14];
   armttbarGL[43]=armttbarGL[12]*armttbarGL[43];
   armttbarGL[42]=3./2.*armttbarGL[43] + 1./4.*armttbarGL[42] + 3*
   armttbarGL[7];
   armttbarGL[40]=11*armttbarGL[11] + armttbarGL[14] + 1./2.*
   armttbarGL[16] + armttbarGL[41] - 145./16. + armttbarGL[40];
   armttbarGL[41]=armttbarGL[3]*armttbarGL[12];
   armttbarGL[40]=3*armttbarGL[40] + 203./4.*armttbarGL[41];
   armttbarGL[40]=armttbarGL[13]*armttbarGL[40];
   armttbarGL[41]=pow(armttbarGL[12],2);
   armttbarGL[43]=armttbarGL[3]*armttbarGL[41];
   armttbarGL[39]=armttbarGL[39] + 1./2.*armttbarGL[40] + 1./2.*
   armttbarGL[42] + armttbarGL[43];
   armttbarGL[38]=1./4.*armttbarGL[39] + armttbarGL[38];
   armttbarGL[38]=MMt*armttbarGL[38];
   armttbarGL[39]= - 51./32.*armttbarGL[12] + armttbarGL[13];
   armttbarGL[39]=armttbarGL[13]*armttbarGL[39];
   armttbarGL[37]=armttbarGL[38] + 1./16.*armttbarGL[37] + 1./128.*
   armttbarGL[41] + armttbarGL[39];

      mttbarGLret = armttbarGL[37]*pow(armttbarGL[4],2)*pow(
      armttbarGL[2],4)*armttbarGL[1];
      return mttbarGLret;
}
