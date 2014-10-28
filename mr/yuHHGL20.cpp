#include <HH.hpp>
std::complex<long double>
HH<OS>::mygl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuHHGL[47], yuHHGLret;

    aryuHHGL[1]=double(boson);
    aryuHHGL[2]=pow(SW,-1);
    aryuHHGL[3]=pow(MMH,-1);
    aryuHHGL[4]=pow(MMW,-1);
    aryuHHGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuHHGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuHHGL[7]=Tsil::I2(0,0,MMH,mu2);
    aryuHHGL[8]=Tsil::I2(0,0,MMt,mu2);
    aryuHHGL[9]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuHHGL[10]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuHHGL[11]=Tsil::A(MMH,mu2);
    aryuHHGL[12]=Tsil::A(MMt,mu2);
    aryuHHGL[13]=std::real(Tsil::B(0,0,MMH,mu2));
    aryuHHGL[14]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuHHGL[15]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuHHGL[16]=Tsil::B(0,MMt,MMH,mu2);
    aryuHHGL[17]=Tsil::B(0,0,MMH,mu2);
    aryuHHGL[18]=Tsil::Beps(MMH,MMH,MMH,mu2);
    aryuHHGL[19]=Tsil::Beps(MMt,MMt,MMH,mu2);
    aryuHHGL[20]=Tsil::Aeps(MMH,mu2);
    aryuHHGL[21]=Tsil::Aeps(MMt,mu2);
    aryuHHGL[22]=prottttt0->Suxv(0);
    aryuHHGL[23]=protHHHHH->M(0);
    aryuHHGL[24]=protHtHtt->M(0);
    aryuHHGL[25]=protttttH->M(0);
    aryuHHGL[26]=protHHHHH->Uzxyv(0);
    aryuHHGL[27]=protHtHtt->Uzxyv(0);
    aryuHHGL[28]=protHtHtt->Uuyxv(0);
    aryuHHGL[29]=protHtHtt->Svyz(0);
    aryuHHGL[30]=prot0H0H0->M(0);
    aryuHHGL[31]=prot0t0tt->M(0);
    aryuHHGL[32]=prot0t0t0->M(0);
    aryuHHGL[33]=prot0000H->M(0);
    aryuHHGL[34]=prot0H0H0->Uuyxv(0);
    aryuHHGL[35]=prot0t0t0->Uuyxv(0);
    aryuHHGL[36]=prot0H0H0->Tyzv(0);
    aryuHHGL[37]=prot0t0t0->Tyzv(0);
    aryuHHGL[38]=1/(4*MMt - MMH);
   aryuHHGL[39]=pow(Pi,2);
   aryuHHGL[39]= - 5 + aryuHHGL[39];
   aryuHHGL[39]=3*aryuHHGL[17] - 9*aryuHHGL[26] - 9*aryuHHGL[34] + 5./2.
   *aryuHHGL[39] + 11*aryuHHGL[36];
   aryuHHGL[40]= - 7 + aryuHHGL[13];
   aryuHHGL[40]=aryuHHGL[13]*aryuHHGL[40];
   aryuHHGL[41]=3*aryuHHGL[13];
   aryuHHGL[42]=21./4.*aryuHHGL[9] + 13./4. + aryuHHGL[41];
   aryuHHGL[42]=aryuHHGL[9]*aryuHHGL[42];
   aryuHHGL[43]= - aryuHHGL[12]*aryuHHGL[38];
   aryuHHGL[44]=aryuHHGL[12]*aryuHHGL[38];
   aryuHHGL[45]=1./8.*aryuHHGL[10] - 1./4. + aryuHHGL[44];
   aryuHHGL[45]=aryuHHGL[10]*aryuHHGL[45];
   aryuHHGL[46]=aryuHHGL[10]*aryuHHGL[38];
   aryuHHGL[46]= - aryuHHGL[38] + 1./2.*aryuHHGL[46];
   aryuHHGL[46]=aryuHHGL[10]*aryuHHGL[46];
   aryuHHGL[46]=1./4.*aryuHHGL[46] + 1./8.*aryuHHGL[38] + 27./2.*
   aryuHHGL[23] + 1./2.*aryuHHGL[33] + 3*aryuHHGL[30];
   aryuHHGL[46]=MMH*aryuHHGL[46];
   aryuHHGL[39]=aryuHHGL[46] + aryuHHGL[45] + aryuHHGL[43] + 3*
   aryuHHGL[42] + 1./4.*aryuHHGL[40] + 1./2.*aryuHHGL[39] + 9*
   aryuHHGL[18];
   aryuHHGL[39]=MMH*aryuHHGL[39];
   aryuHHGL[40]= - 9*aryuHHGL[9];
   aryuHHGL[42]=aryuHHGL[40] - 5 + aryuHHGL[41];
   aryuHHGL[42]=1./2.*aryuHHGL[42] + aryuHHGL[44];
   aryuHHGL[42]=aryuHHGL[12]*aryuHHGL[42];
   aryuHHGL[43]=aryuHHGL[7] + 3*aryuHHGL[5];
   aryuHHGL[43]=1./8.*aryuHHGL[20] + 1./8.*aryuHHGL[43] + aryuHHGL[21];
   aryuHHGL[44]=45./16.*aryuHHGL[9] - 1 + 13./16.*aryuHHGL[13];
   aryuHHGL[44]=aryuHHGL[11]*aryuHHGL[44];
   aryuHHGL[45]= - aryuHHGL[10]*aryuHHGL[12];
   aryuHHGL[39]=1./8.*aryuHHGL[39] + 1./8.*aryuHHGL[45] + 1./4.*
   aryuHHGL[42] + 1./2.*aryuHHGL[43] + aryuHHGL[44];
   aryuHHGL[39]=MMH*aryuHHGL[39];
   aryuHHGL[42]= - 3*aryuHHGL[11] - aryuHHGL[12];
   aryuHHGL[42]=aryuHHGL[12]*aryuHHGL[42];
   aryuHHGL[43]=pow(aryuHHGL[11],2);
   aryuHHGL[42]=51./8.*aryuHHGL[43] + aryuHHGL[42];
   aryuHHGL[39]=1./4.*aryuHHGL[42] + aryuHHGL[39];
   aryuHHGL[42]=7./4. - aryuHHGL[28];
   aryuHHGL[42]=1./8.*aryuHHGL[17] + 1./4.*aryuHHGL[42] - aryuHHGL[27];
   aryuHHGL[41]= - aryuHHGL[14] + aryuHHGL[41];
   aryuHHGL[41]= - 1./16.*aryuHHGL[10] + 1./4.*aryuHHGL[41] + 3*
   aryuHHGL[9];
   aryuHHGL[41]=aryuHHGL[10]*aryuHHGL[41];
   aryuHHGL[41]=aryuHHGL[41] + 27./16.*aryuHHGL[9] - 1./16.*
   aryuHHGL[13] - 3./4.*aryuHHGL[16] + 3*aryuHHGL[18] + 3./4.*
   aryuHHGL[19] + 3*aryuHHGL[42] + 1./2.*aryuHHGL[37];
   aryuHHGL[41]=MMH*aryuHHGL[41];
   aryuHHGL[42]=5./2. + aryuHHGL[14];
   aryuHHGL[42]= - 9./2.*aryuHHGL[9] + 5./2.*aryuHHGL[42] - 
   aryuHHGL[13];
   aryuHHGL[42]=aryuHHGL[12]*aryuHHGL[42];
   aryuHHGL[43]=5*aryuHHGL[11] + 3*aryuHHGL[12];
   aryuHHGL[43]=aryuHHGL[10]*aryuHHGL[43];
   aryuHHGL[44]= - 11*aryuHHGL[11] + 13*aryuHHGL[12];
   aryuHHGL[44]=aryuHHGL[3]*aryuHHGL[12]*aryuHHGL[44];
   aryuHHGL[41]=aryuHHGL[41] + 1./4.*aryuHHGL[44] + 1./4.*aryuHHGL[43]
    + aryuHHGL[42] + 3./8.*aryuHHGL[11] + 3./2.*aryuHHGL[20] + 5*
   aryuHHGL[21] - aryuHHGL[22] - 3./2.*aryuHHGL[6];
   aryuHHGL[42]=7*aryuHHGL[8] - 17*aryuHHGL[29] - aryuHHGL[22];
   aryuHHGL[42]= - 15./2.*aryuHHGL[11] - 5./2.*aryuHHGL[20] - 9*
   aryuHHGL[21] + 1./2.*aryuHHGL[42] + 11*aryuHHGL[6];
   aryuHHGL[43]=21./4. - aryuHHGL[15];
   aryuHHGL[43]=1./4.*aryuHHGL[43] - aryuHHGL[14];
   aryuHHGL[43]=aryuHHGL[12]*aryuHHGL[43];
   aryuHHGL[44]= - 2*aryuHHGL[11];
   aryuHHGL[45]=aryuHHGL[44] + 1./4.*aryuHHGL[12];
   aryuHHGL[45]=aryuHHGL[10]*aryuHHGL[45];
   aryuHHGL[44]=aryuHHGL[44] - aryuHHGL[12];
   aryuHHGL[44]=aryuHHGL[3]*aryuHHGL[12]*aryuHHGL[44];
   aryuHHGL[42]=2*aryuHHGL[44] + aryuHHGL[45] + 1./4.*aryuHHGL[42] + 5*
   aryuHHGL[43];
   aryuHHGL[42]=aryuHHGL[3]*aryuHHGL[42];
   aryuHHGL[43]=3./8. - aryuHHGL[35];
   aryuHHGL[43]=aryuHHGL[40] + aryuHHGL[14] + 9./4.*aryuHHGL[16] - 9*
   aryuHHGL[18] + 5./2.*aryuHHGL[19] + 5./2.*aryuHHGL[43] + 9*
   aryuHHGL[27];
   aryuHHGL[44]= - 5 + aryuHHGL[15];
   aryuHHGL[44]=1./2.*aryuHHGL[44] + 7*aryuHHGL[14];
   aryuHHGL[40]= - 1./8.*aryuHHGL[10] + aryuHHGL[40] + 1./2.*
   aryuHHGL[44] - aryuHHGL[13];
   aryuHHGL[40]=aryuHHGL[10]*aryuHHGL[40];
   aryuHHGL[40]=1./2.*aryuHHGL[43] + aryuHHGL[40];
   aryuHHGL[43]= - 1 - 5*aryuHHGL[15];
   aryuHHGL[43]= - 7./8.*aryuHHGL[10] + 1./4.*aryuHHGL[43] - 5*
   aryuHHGL[14];
   aryuHHGL[43]=aryuHHGL[10]*aryuHHGL[43];
   aryuHHGL[44]= - aryuHHGL[10]*aryuHHGL[11];
   aryuHHGL[44]=2*aryuHHGL[44] + 9./2.*aryuHHGL[12] - 4*aryuHHGL[11] + 
   2*aryuHHGL[20] + 4*aryuHHGL[21] - 1./2.*aryuHHGL[8] - 2*aryuHHGL[6];
   aryuHHGL[44]=aryuHHGL[3]*aryuHHGL[44];
   aryuHHGL[43]=aryuHHGL[44] + aryuHHGL[43] - aryuHHGL[14] - 1./4.*
   aryuHHGL[15] + 5./16.*aryuHHGL[16] - 13./4.*aryuHHGL[19] - 3./8.*
   aryuHHGL[37] + 2*aryuHHGL[28] + 6 + 5./4.*aryuHHGL[35];
   aryuHHGL[43]=aryuHHGL[3]*aryuHHGL[43];
   aryuHHGL[44]= - aryuHHGL[16] - 9 - aryuHHGL[37];
   aryuHHGL[44]=aryuHHGL[3]*aryuHHGL[44];
   aryuHHGL[44]= - 2*aryuHHGL[25] + 1./2.*aryuHHGL[44];
   aryuHHGL[44]=MMt*aryuHHGL[3]*aryuHHGL[44];
   aryuHHGL[43]=aryuHHGL[44] + aryuHHGL[43] + aryuHHGL[25] - 1./2.*
   aryuHHGL[32] - 6*aryuHHGL[24];
   aryuHHGL[43]=MMt*aryuHHGL[43];
   aryuHHGL[44]=1./2.*aryuHHGL[25] - aryuHHGL[31] + 9*aryuHHGL[24];
   aryuHHGL[44]=MMH*aryuHHGL[44];
   aryuHHGL[40]=aryuHHGL[43] + 1./4.*aryuHHGL[44] + 1./2.*aryuHHGL[40]
    + aryuHHGL[42];
   aryuHHGL[40]=MMt*aryuHHGL[40];
   aryuHHGL[40]=1./2.*aryuHHGL[41] + aryuHHGL[40];
   aryuHHGL[40]=MMt*aryuHHGL[40];
   aryuHHGL[39]=1./2.*aryuHHGL[39] + aryuHHGL[40];

      yuHHGLret = 3*aryuHHGL[39]*pow(aryuHHGL[4],2)*pow(aryuHHGL[2],4)*
      aryuHHGL[1];
      return yuHHGLret;
}
