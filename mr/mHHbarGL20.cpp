#include <HH.hpp>
std::complex<long double>
HH<OS>::mgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbarGL[54], mHHbarGLret;

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
   armHHbarGL[39]=armHHbarGL[16] + armHHbarGL[28];
   armHHbarGL[40]=armHHbarGL[18] - armHHbarGL[27];
   armHHbarGL[41]=MMt*armHHbarGL[24];
   armHHbarGL[42]=1./4.*armHHbarGL[10];
   armHHbarGL[43]= - armHHbarGL[14] - armHHbarGL[42];
   armHHbarGL[43]=armHHbarGL[10]*armHHbarGL[43];
   armHHbarGL[43]=21./4. + armHHbarGL[43];
   armHHbarGL[39]=9./2.*armHHbarGL[41] + 3./8.*armHHbarGL[17] + 1./4.*
   armHHbarGL[43] + 3*armHHbarGL[40] - 3./4.*armHHbarGL[39];
   armHHbarGL[39]=MMt*armHHbarGL[39];
   armHHbarGL[40]=1./4.*armHHbarGL[13];
   armHHbarGL[43]= - 1./4. + 3*armHHbarGL[10];
   armHHbarGL[43]=MMt*armHHbarGL[43];
   armHHbarGL[43]=armHHbarGL[43] + 13./4.*armHHbarGL[11];
   armHHbarGL[43]=armHHbarGL[43]*armHHbarGL[40];
   armHHbarGL[44]=3*armHHbarGL[9];
   armHHbarGL[45]=9./16. + armHHbarGL[10];
   armHHbarGL[45]=MMt*armHHbarGL[45];
   armHHbarGL[45]=armHHbarGL[45] + 15./16.*armHHbarGL[11];
   armHHbarGL[45]=armHHbarGL[45]*armHHbarGL[44];
   armHHbarGL[46]=armHHbarGL[20] + armHHbarGL[7];
   armHHbarGL[47]=pow(MMt,2);
   armHHbarGL[48]=armHHbarGL[25]*armHHbarGL[47];
   armHHbarGL[49]=3*armHHbarGL[13];
   armHHbarGL[50]= - 9*armHHbarGL[9] + armHHbarGL[49] - 5 - 
   armHHbarGL[10];
   armHHbarGL[50]=armHHbarGL[12]*armHHbarGL[50];
   armHHbarGL[51]=3./4.*MMt;
   armHHbarGL[52]=armHHbarGL[51]*armHHbarGL[19];
   armHHbarGL[39]=armHHbarGL[52] + 1./8.*armHHbarGL[50] + 
   armHHbarGL[45] + armHHbarGL[43] + 1./4.*armHHbarGL[48] + 1./2.*
   armHHbarGL[21] - armHHbarGL[11] + 3./16.*armHHbarGL[5] + 
   armHHbarGL[39] + 1./16.*armHHbarGL[46];
   armHHbarGL[39]=armHHbarGL[1]*armHHbarGL[39];
   armHHbarGL[43]=armHHbarGL[26] + armHHbarGL[34];
   armHHbarGL[45]=1./2.*armHHbarGL[10];
   armHHbarGL[46]=armHHbarGL[45] - 1;
   armHHbarGL[46]=armHHbarGL[46]*armHHbarGL[10];
   armHHbarGL[48]= - 25 + armHHbarGL[46];
   armHHbarGL[43]=1./2.*armHHbarGL[48] + 3*armHHbarGL[17] - 9*
   armHHbarGL[43];
   armHHbarGL[48]= - 7 + armHHbarGL[13];
   armHHbarGL[40]=armHHbarGL[48]*armHHbarGL[40];
   armHHbarGL[48]=21./4.*armHHbarGL[9] + 13./4. + armHHbarGL[49];
   armHHbarGL[44]=armHHbarGL[48]*armHHbarGL[44];
   armHHbarGL[48]=pow(Pi,2);
   armHHbarGL[49]= - 1 + armHHbarGL[10];
   armHHbarGL[49]=armHHbarGL[12]*armHHbarGL[49]*armHHbarGL[38];
   armHHbarGL[40]=armHHbarGL[49] + armHHbarGL[44] + armHHbarGL[40] + 11.
   /2.*armHHbarGL[36] + 5./4.*armHHbarGL[48] + 1./2.*armHHbarGL[43] + 9
   *armHHbarGL[18];
   armHHbarGL[40]=armHHbarGL[1]*armHHbarGL[40];
   armHHbarGL[43]=1./4.*armHHbarGL[1];
   armHHbarGL[44]=armHHbarGL[43]*armHHbarGL[38];
   armHHbarGL[46]=1./2. + armHHbarGL[46];
   armHHbarGL[46]=armHHbarGL[46]*armHHbarGL[44];
   armHHbarGL[48]=27./2.*armHHbarGL[23] + 3*armHHbarGL[30] + 1./2.*
   armHHbarGL[33];
   armHHbarGL[48]=armHHbarGL[1]*armHHbarGL[48];
   armHHbarGL[46]=armHHbarGL[46] + armHHbarGL[48];
   armHHbarGL[46]=MMH*armHHbarGL[46];
   armHHbarGL[40]=armHHbarGL[46] + armHHbarGL[40];
   armHHbarGL[40]=MMH*armHHbarGL[40];
   armHHbarGL[46]=1./2.*armHHbarGL[1];
   armHHbarGL[48]= - armHHbarGL[31]*armHHbarGL[47]*armHHbarGL[46];
   armHHbarGL[44]=pow(armHHbarGL[12],2)*armHHbarGL[44];
   armHHbarGL[39]=1./8.*armHHbarGL[40] + armHHbarGL[48] + 
   armHHbarGL[44] + armHHbarGL[39];
   armHHbarGL[39]=MMH*armHHbarGL[39];
   armHHbarGL[40]=pow(MMt,3);
   armHHbarGL[44]=pow(MMt,4);
   armHHbarGL[48]= - armHHbarGL[3]*armHHbarGL[44];
   armHHbarGL[48]= - 3./4.*armHHbarGL[40] + armHHbarGL[48];
   armHHbarGL[48]=armHHbarGL[1]*armHHbarGL[3]*armHHbarGL[48];
   armHHbarGL[49]=MMH*MMt*armHHbarGL[46];
   armHHbarGL[48]=armHHbarGL[48] + armHHbarGL[49];
   armHHbarGL[48]=armHHbarGL[37]*armHHbarGL[48];
   armHHbarGL[39]=armHHbarGL[39] + armHHbarGL[48];
   armHHbarGL[48]=armHHbarGL[20] - armHHbarGL[6];
   armHHbarGL[49]= - 2 - armHHbarGL[10];
   armHHbarGL[49]=armHHbarGL[11]*armHHbarGL[49];
   armHHbarGL[49]=armHHbarGL[48] + armHHbarGL[49];
   armHHbarGL[49]= - 1./2.*armHHbarGL[8] + 4*armHHbarGL[21] + 2*
   armHHbarGL[49];
   armHHbarGL[49]=armHHbarGL[40]*armHHbarGL[49];
   armHHbarGL[50]= - 9 - armHHbarGL[16];
   armHHbarGL[50]=armHHbarGL[50]*armHHbarGL[44];
   armHHbarGL[52]= - 2*armHHbarGL[12] - 4*armHHbarGL[11];
   armHHbarGL[52]=armHHbarGL[47]*armHHbarGL[52];
   armHHbarGL[52]=9./2.*armHHbarGL[40] + armHHbarGL[52];
   armHHbarGL[52]=armHHbarGL[12]*armHHbarGL[52];
   armHHbarGL[49]=armHHbarGL[52] + 1./2.*armHHbarGL[50] + 
   armHHbarGL[49];
   armHHbarGL[49]=armHHbarGL[3]*armHHbarGL[49];
   armHHbarGL[50]= - 7./8.*armHHbarGL[10] - 1./4. - 5*armHHbarGL[14];
   armHHbarGL[50]=armHHbarGL[10]*armHHbarGL[50];
   armHHbarGL[50]=5./16.*armHHbarGL[16] + 2*armHHbarGL[28] + 
   armHHbarGL[50] + 6 - armHHbarGL[14];
   armHHbarGL[50]=MMt*armHHbarGL[50];
   armHHbarGL[52]= - 15./8. - 2*armHHbarGL[10];
   armHHbarGL[52]=armHHbarGL[11]*armHHbarGL[52];
   armHHbarGL[50]= - 5./8.*armHHbarGL[20] + 7./8.*armHHbarGL[8] - 1./8.
   *armHHbarGL[22] - 9./4.*armHHbarGL[21] + armHHbarGL[52] + 11./4.*
   armHHbarGL[6] - 17./8.*armHHbarGL[29] + armHHbarGL[50];
   armHHbarGL[50]=armHHbarGL[47]*armHHbarGL[50];
   armHHbarGL[52]=21./16. - armHHbarGL[14];
   armHHbarGL[52]=5*armHHbarGL[52] + armHHbarGL[42];
   armHHbarGL[52]=armHHbarGL[52]*armHHbarGL[47];
   armHHbarGL[53]=13./8.*armHHbarGL[12] - 11./8.*armHHbarGL[11];
   armHHbarGL[53]=MMt*armHHbarGL[53];
   armHHbarGL[52]=armHHbarGL[52] + armHHbarGL[53];
   armHHbarGL[52]=armHHbarGL[12]*armHHbarGL[52];
   armHHbarGL[44]=armHHbarGL[25]*armHHbarGL[44];
   armHHbarGL[53]=armHHbarGL[35]*armHHbarGL[40];
   armHHbarGL[44]=armHHbarGL[49] + armHHbarGL[52] + 5./4.*
   armHHbarGL[53] - 2*armHHbarGL[44] + armHHbarGL[50];
   armHHbarGL[44]=armHHbarGL[3]*armHHbarGL[44];
   armHHbarGL[49]=9./2.*armHHbarGL[9];
   armHHbarGL[50]= - 1./2. - armHHbarGL[10];
   armHHbarGL[50]=armHHbarGL[50]*armHHbarGL[49];
   armHHbarGL[42]= - armHHbarGL[42] - 5./2. + 7*armHHbarGL[14];
   armHHbarGL[42]=armHHbarGL[10]*armHHbarGL[42];
   armHHbarGL[42]=9./4.*armHHbarGL[16] + 9*armHHbarGL[27] + 
   armHHbarGL[42] + 15./16. + armHHbarGL[14];
   armHHbarGL[52]= - armHHbarGL[13]*armHHbarGL[45];
   armHHbarGL[41]=armHHbarGL[50] + armHHbarGL[52] - 9./4.*
   armHHbarGL[18] - 5./8.*armHHbarGL[35] + 1./4.*armHHbarGL[42] - 6*
   armHHbarGL[41];
   armHHbarGL[41]=armHHbarGL[47]*armHHbarGL[41];
   armHHbarGL[42]=armHHbarGL[51]*armHHbarGL[48];
   armHHbarGL[48]= - armHHbarGL[49] - armHHbarGL[13];
   armHHbarGL[48]=MMt*armHHbarGL[48];
   armHHbarGL[49]=5./2. + armHHbarGL[14];
   armHHbarGL[49]=5*armHHbarGL[49] + 3./2.*armHHbarGL[10];
   armHHbarGL[49]=MMt*armHHbarGL[49];
   armHHbarGL[49]=armHHbarGL[49] - 3./2.*armHHbarGL[11];
   armHHbarGL[48]= - 1./4.*armHHbarGL[12] + 1./2.*armHHbarGL[49] + 
   armHHbarGL[48];
   armHHbarGL[48]=armHHbarGL[12]*armHHbarGL[48];
   armHHbarGL[49]= - 1./2.*armHHbarGL[22] + 5./2.*armHHbarGL[21];
   armHHbarGL[49]=MMt*armHHbarGL[49];
   armHHbarGL[50]=5*armHHbarGL[10];
   armHHbarGL[51]=3./2. + armHHbarGL[50];
   armHHbarGL[51]=MMt*armHHbarGL[51];
   armHHbarGL[51]=armHHbarGL[51] + 51./8.*armHHbarGL[11];
   armHHbarGL[51]=armHHbarGL[11]*armHHbarGL[51];
   armHHbarGL[52]=armHHbarGL[25]*armHHbarGL[40];
   armHHbarGL[41]=armHHbarGL[44] + 1./2.*armHHbarGL[48] + 
   armHHbarGL[52] + 1./8.*armHHbarGL[51] + armHHbarGL[41] + 
   armHHbarGL[42] + armHHbarGL[49];
   armHHbarGL[41]=armHHbarGL[1]*armHHbarGL[41];
   armHHbarGL[42]=armHHbarGL[47]*armHHbarGL[45];
   armHHbarGL[44]= - 1 - armHHbarGL[50];
   armHHbarGL[44]=armHHbarGL[44]*armHHbarGL[40];
   armHHbarGL[45]=armHHbarGL[12]*armHHbarGL[47];
   armHHbarGL[44]=armHHbarGL[44] - 5*armHHbarGL[45];
   armHHbarGL[44]=armHHbarGL[3]*armHHbarGL[44];
   armHHbarGL[42]=armHHbarGL[42] + armHHbarGL[44];
   armHHbarGL[42]=armHHbarGL[15]*armHHbarGL[42];
   armHHbarGL[44]=armHHbarGL[3]*armHHbarGL[40];
   armHHbarGL[44]=5./2.*armHHbarGL[47] - 13*armHHbarGL[44];
   armHHbarGL[44]=armHHbarGL[19]*armHHbarGL[44];
   armHHbarGL[42]=armHHbarGL[44] + armHHbarGL[42];
   armHHbarGL[42]=armHHbarGL[43]*armHHbarGL[42];
   armHHbarGL[40]= - armHHbarGL[32]*armHHbarGL[40]*armHHbarGL[46];
   armHHbarGL[39]=armHHbarGL[41] + armHHbarGL[40] + armHHbarGL[42] + 1./
   2.*armHHbarGL[39];

      mHHbarGLret = 3*armHHbarGL[39]*pow(armHHbarGL[4],2)*pow(
      armHHbarGL[2],4);
      return mHHbarGLret;
}
