#include <HH.hpp>
std::complex<long double>
HH<OS>::mgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbarGL[47], mHHbarGLret;

    armHHbarGL[1]=double(boson);
    armHHbarGL[2]=pow(SW,-1);
    armHHbarGL[3]=pow(MMH,-1);
    armHHbarGL[4]=pow(MMW,-1);
    armHHbarGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    armHHbarGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    armHHbarGL[7]=Tsil::I2(0,0,MMH,mu2);
    armHHbarGL[8]=Tsil::I2(0,0,MMt,mu2);
    armHHbarGL[9]=Tsil::B(MMH,MMH,MMH,mu2);
    armHHbarGL[10]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbarGL[11]=Tsil::A(MMH,mu2);
    armHHbarGL[12]=Tsil::A(MMt,mu2);
    armHHbarGL[13]=std::real(Tsil::B(0,0,MMH,mu2));
    armHHbarGL[14]=Tsil::B(MMH,MMt,MMt,mu2);
    armHHbarGL[15]=std::real(Tsil::B(0,0,MMt,mu2));
    armHHbarGL[16]=Tsil::B(0,MMt,MMH,mu2);
    armHHbarGL[17]=Tsil::B(0,0,MMH,mu2);
    armHHbarGL[18]=Tsil::Beps(MMH,MMH,MMH,mu2);
    armHHbarGL[19]=Tsil::Beps(MMt,MMt,MMH,mu2);
    armHHbarGL[20]=Tsil::Aeps(MMH,mu2);
    armHHbarGL[21]=Tsil::Aeps(MMt,mu2);
    armHHbarGL[22]=prottttt0->Suxv(0);
    armHHbarGL[23]=protHHHHH->M(0);
    armHHbarGL[24]=protHtHtt->M(0);
    armHHbarGL[25]=protttttH->M(0);
    armHHbarGL[26]=protHHHHH->Uzxyv(0);
    armHHbarGL[27]=protHtHtt->Uzxyv(0);
    armHHbarGL[28]=protHtHtt->Uuyxv(0);
    armHHbarGL[29]=protHtHtt->Svyz(0);
    armHHbarGL[30]=prot0H0H0->M(0);
    armHHbarGL[31]=prot0t0tt->M(0);
    armHHbarGL[32]=prot0t0t0->M(0);
    armHHbarGL[33]=prot0000H->M(0);
    armHHbarGL[34]=prot0H0H0->Uuyxv(0);
    armHHbarGL[35]=prot0t0t0->Uuyxv(0);
    armHHbarGL[36]=prot0H0H0->Tyzv(0);
    armHHbarGL[37]=prot0t0t0->Tyzv(0);
    armHHbarGL[38]=1/(4*MMt - MMH);
   armHHbarGL[39]=pow(Pi,2);
   armHHbarGL[39]= - 5 + armHHbarGL[39];
   armHHbarGL[39]=3*armHHbarGL[17] - 9*armHHbarGL[26] - 9*
   armHHbarGL[34] + 5./2.*armHHbarGL[39] + 11*armHHbarGL[36];
   armHHbarGL[40]= - 7 + armHHbarGL[13];
   armHHbarGL[40]=armHHbarGL[13]*armHHbarGL[40];
   armHHbarGL[41]=3*armHHbarGL[13];
   armHHbarGL[42]=21./4.*armHHbarGL[9] + 13./4. + armHHbarGL[41];
   armHHbarGL[42]=armHHbarGL[9]*armHHbarGL[42];
   armHHbarGL[43]= - armHHbarGL[12]*armHHbarGL[38];
   armHHbarGL[44]=armHHbarGL[12]*armHHbarGL[38];
   armHHbarGL[45]=1./8.*armHHbarGL[10] - 1./4. + armHHbarGL[44];
   armHHbarGL[45]=armHHbarGL[10]*armHHbarGL[45];
   armHHbarGL[46]=armHHbarGL[10]*armHHbarGL[38];
   armHHbarGL[46]= - armHHbarGL[38] + 1./2.*armHHbarGL[46];
   armHHbarGL[46]=armHHbarGL[10]*armHHbarGL[46];
   armHHbarGL[46]=1./4.*armHHbarGL[46] + 1./8.*armHHbarGL[38] + 27./2.*
   armHHbarGL[23] + 1./2.*armHHbarGL[33] + 3*armHHbarGL[30];
   armHHbarGL[46]=MMH*armHHbarGL[46];
   armHHbarGL[39]=armHHbarGL[46] + armHHbarGL[45] + armHHbarGL[43] + 3*
   armHHbarGL[42] + 1./4.*armHHbarGL[40] + 1./2.*armHHbarGL[39] + 9*
   armHHbarGL[18];
   armHHbarGL[39]=MMH*armHHbarGL[39];
   armHHbarGL[40]= - 9*armHHbarGL[9];
   armHHbarGL[42]=armHHbarGL[40] - 5 + armHHbarGL[41];
   armHHbarGL[42]=1./2.*armHHbarGL[42] + armHHbarGL[44];
   armHHbarGL[42]=armHHbarGL[12]*armHHbarGL[42];
   armHHbarGL[43]=armHHbarGL[7] + 3*armHHbarGL[5];
   armHHbarGL[43]=1./8.*armHHbarGL[20] + 1./8.*armHHbarGL[43] + 
   armHHbarGL[21];
   armHHbarGL[44]=45./16.*armHHbarGL[9] - 1 + 13./16.*armHHbarGL[13];
   armHHbarGL[44]=armHHbarGL[11]*armHHbarGL[44];
   armHHbarGL[45]= - armHHbarGL[10]*armHHbarGL[12];
   armHHbarGL[39]=1./8.*armHHbarGL[39] + 1./8.*armHHbarGL[45] + 1./4.*
   armHHbarGL[42] + 1./2.*armHHbarGL[43] + armHHbarGL[44];
   armHHbarGL[39]=MMH*armHHbarGL[39];
   armHHbarGL[42]= - 3*armHHbarGL[11] - armHHbarGL[12];
   armHHbarGL[42]=armHHbarGL[12]*armHHbarGL[42];
   armHHbarGL[43]=pow(armHHbarGL[11],2);
   armHHbarGL[42]=51./8.*armHHbarGL[43] + armHHbarGL[42];
   armHHbarGL[39]=1./4.*armHHbarGL[42] + armHHbarGL[39];
   armHHbarGL[42]=7./4. - armHHbarGL[28];
   armHHbarGL[42]=1./8.*armHHbarGL[17] + 1./4.*armHHbarGL[42] - 
   armHHbarGL[27];
   armHHbarGL[41]= - armHHbarGL[14] + armHHbarGL[41];
   armHHbarGL[41]= - 1./16.*armHHbarGL[10] + 1./4.*armHHbarGL[41] + 3*
   armHHbarGL[9];
   armHHbarGL[41]=armHHbarGL[10]*armHHbarGL[41];
   armHHbarGL[41]=armHHbarGL[41] + 27./16.*armHHbarGL[9] - 1./16.*
   armHHbarGL[13] - 3./4.*armHHbarGL[16] + 3*armHHbarGL[18] + 3./4.*
   armHHbarGL[19] + 3*armHHbarGL[42] + 1./2.*armHHbarGL[37];
   armHHbarGL[41]=MMH*armHHbarGL[41];
   armHHbarGL[42]=5./2. + armHHbarGL[14];
   armHHbarGL[42]= - 9./2.*armHHbarGL[9] + 5./2.*armHHbarGL[42] - 
   armHHbarGL[13];
   armHHbarGL[42]=armHHbarGL[12]*armHHbarGL[42];
   armHHbarGL[43]=5*armHHbarGL[11] + 3*armHHbarGL[12];
   armHHbarGL[43]=armHHbarGL[10]*armHHbarGL[43];
   armHHbarGL[44]= - 11*armHHbarGL[11] + 13*armHHbarGL[12];
   armHHbarGL[44]=armHHbarGL[3]*armHHbarGL[12]*armHHbarGL[44];
   armHHbarGL[41]=armHHbarGL[41] + 1./4.*armHHbarGL[44] + 1./4.*
   armHHbarGL[43] + armHHbarGL[42] + 3./8.*armHHbarGL[11] + 3./2.*
   armHHbarGL[20] + 5*armHHbarGL[21] - armHHbarGL[22] - 3./2.*
   armHHbarGL[6];
   armHHbarGL[42]=7*armHHbarGL[8] - 17*armHHbarGL[29] - armHHbarGL[22];
   armHHbarGL[42]= - 15./2.*armHHbarGL[11] - 5./2.*armHHbarGL[20] - 9*
   armHHbarGL[21] + 1./2.*armHHbarGL[42] + 11*armHHbarGL[6];
   armHHbarGL[43]=21./4. - armHHbarGL[15];
   armHHbarGL[43]=1./4.*armHHbarGL[43] - armHHbarGL[14];
   armHHbarGL[43]=armHHbarGL[12]*armHHbarGL[43];
   armHHbarGL[44]= - 2*armHHbarGL[11];
   armHHbarGL[45]=armHHbarGL[44] + 1./4.*armHHbarGL[12];
   armHHbarGL[45]=armHHbarGL[10]*armHHbarGL[45];
   armHHbarGL[44]=armHHbarGL[44] - armHHbarGL[12];
   armHHbarGL[44]=armHHbarGL[3]*armHHbarGL[12]*armHHbarGL[44];
   armHHbarGL[42]=2*armHHbarGL[44] + armHHbarGL[45] + 1./4.*
   armHHbarGL[42] + 5*armHHbarGL[43];
   armHHbarGL[42]=armHHbarGL[3]*armHHbarGL[42];
   armHHbarGL[43]=3./8. - armHHbarGL[35];
   armHHbarGL[43]=armHHbarGL[40] + armHHbarGL[14] + 9./4.*
   armHHbarGL[16] - 9*armHHbarGL[18] + 5./2.*armHHbarGL[19] + 5./2.*
   armHHbarGL[43] + 9*armHHbarGL[27];
   armHHbarGL[44]= - 5 + armHHbarGL[15];
   armHHbarGL[44]=1./2.*armHHbarGL[44] + 7*armHHbarGL[14];
   armHHbarGL[40]= - 1./8.*armHHbarGL[10] + armHHbarGL[40] + 1./2.*
   armHHbarGL[44] - armHHbarGL[13];
   armHHbarGL[40]=armHHbarGL[10]*armHHbarGL[40];
   armHHbarGL[40]=1./2.*armHHbarGL[43] + armHHbarGL[40];
   armHHbarGL[43]= - 1 - 5*armHHbarGL[15];
   armHHbarGL[43]= - 7./8.*armHHbarGL[10] + 1./4.*armHHbarGL[43] - 5*
   armHHbarGL[14];
   armHHbarGL[43]=armHHbarGL[10]*armHHbarGL[43];
   armHHbarGL[44]= - armHHbarGL[10]*armHHbarGL[11];
   armHHbarGL[44]=2*armHHbarGL[44] + 9./2.*armHHbarGL[12] - 4*
   armHHbarGL[11] + 2*armHHbarGL[20] + 4*armHHbarGL[21] - 1./2.*
   armHHbarGL[8] - 2*armHHbarGL[6];
   armHHbarGL[44]=armHHbarGL[3]*armHHbarGL[44];
   armHHbarGL[43]=armHHbarGL[44] + armHHbarGL[43] - armHHbarGL[14] - 1./
   4.*armHHbarGL[15] + 5./16.*armHHbarGL[16] - 13./4.*armHHbarGL[19] - 
   3./8.*armHHbarGL[37] + 2*armHHbarGL[28] + 6 + 5./4.*armHHbarGL[35];
   armHHbarGL[43]=armHHbarGL[3]*armHHbarGL[43];
   armHHbarGL[44]= - armHHbarGL[16] - 9 - armHHbarGL[37];
   armHHbarGL[44]=armHHbarGL[3]*armHHbarGL[44];
   armHHbarGL[44]= - 2*armHHbarGL[25] + 1./2.*armHHbarGL[44];
   armHHbarGL[44]=MMt*armHHbarGL[3]*armHHbarGL[44];
   armHHbarGL[43]=armHHbarGL[44] + armHHbarGL[43] + armHHbarGL[25] - 1./
   2.*armHHbarGL[32] - 6*armHHbarGL[24];
   armHHbarGL[43]=MMt*armHHbarGL[43];
   armHHbarGL[44]=1./2.*armHHbarGL[25] - armHHbarGL[31] + 9*
   armHHbarGL[24];
   armHHbarGL[44]=MMH*armHHbarGL[44];
   armHHbarGL[40]=armHHbarGL[43] + 1./4.*armHHbarGL[44] + 1./2.*
   armHHbarGL[40] + armHHbarGL[42];
   armHHbarGL[40]=MMt*armHHbarGL[40];
   armHHbarGL[40]=1./2.*armHHbarGL[41] + armHHbarGL[40];
   armHHbarGL[40]=MMt*armHHbarGL[40];
   armHHbarGL[39]=1./2.*armHHbarGL[39] + armHHbarGL[40];

      mHHbarGLret = 3*armHHbarGL[39]*pow(armHHbarGL[4],2)*pow(
      armHHbarGL[2],4)*armHHbarGL[1];
      return mHHbarGLret;
}
