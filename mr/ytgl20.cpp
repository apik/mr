#include <tt.hpp>
std::complex<long double> tt::mygl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[55];

    mytgl[1]=pow(SW,-1);
    mytgl[2]=pow(MMH,-1);
    mytgl[3]=pow(MMW,-1);
    mytgl[4]=pow(MMt,-1);
    mytgl[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    mytgl[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    mytgl[7]=Tsil::I2(0,MMH,MMt,mu2);
    mytgl[8]=Tsil::I2(0,0,MMH,mu2);
    mytgl[9]=Tsil::I2(0,0,MMt,mu2);
    mytgl[10]=Tsil::B(MMH,MMH,MMH,mu2);
    mytgl[11]=Tsil::B(MMH,MMt,MMt,mu2);
    mytgl[12]=Tsil::A(MMH,mu2);
    mytgl[13]=Tsil::A(MMt,mu2);
    mytgl[14]=Tsil::B(MMt,MMt,MMH,mu2);
    mytgl[15]=std::real(Tsil::B(0,0,MMH,mu2));
    mytgl[16]=std::real(Tsil::B(0,0,MMt,mu2));
    mytgl[17]=Tsil::B(0,MMH,MMt,mu2);
    mytgl[18]=Tsil::Beps(MMH,MMt,MMt,mu2);
    mytgl[19]=Tsil::Aeps(MMH,mu2);
    mytgl[20]=Tsil::Aeps(MMt,mu2);
    mytgl[21]=prot0ttHt->Tuxv(0);
    mytgl[22]=protHHttH->M(0);
    mytgl[23]=protHttHt->M(0);
    mytgl[24]=protHHttH->Uxzuv(0);
    mytgl[25]=protHttHt->Txuv(0);
    mytgl[26]=protHHttH->Suxv(0);
    mytgl[27]=protHHttH->Uuyxv(0);
    mytgl[28]=protWt000->Tyzv(0);
    mytgl[29]=protH0tt0->M(0);
    mytgl[30]=protH0t00->M(0);
    mytgl[31]=prot0Htt0->M(0);
    mytgl[32]=prot0H0t0->M(0);
    mytgl[33]=prot00ttH->M(0);
    mytgl[34]=protH0t00->Uxzuv(0);
    mytgl[35]=protH0t00->Uzxyv(0);
    mytgl[36]=protH0tt0->Uuyxv(0);
    mytgl[37]=protH0t00->Txuv(0);
   mytgl[38]=1./2.*mytgl[11];
   mytgl[39]= - 61./4. + mytgl[11];
   mytgl[39]=mytgl[39]*mytgl[38];
   mytgl[40]= - mytgl[37] - 1 - mytgl[17];
   mytgl[41]=1./2.*MMH;
   mytgl[40]=mytgl[41]*mytgl[4]*mytgl[40];
   mytgl[42]=mytgl[25] - mytgl[27];
   mytgl[43]=mytgl[13] + mytgl[20];
   mytgl[43]=mytgl[4]*mytgl[43];
   mytgl[43]=mytgl[43] + mytgl[36];
   mytgl[44]=mytgl[11] - 1./4.;
   mytgl[45]=mytgl[15]*mytgl[44];
   mytgl[45]=mytgl[45] - mytgl[35];
   mytgl[46]=pow(Pi,2);
   mytgl[39]=mytgl[40] - 23./4.*mytgl[21] - 21./4.*mytgl[18] + 3./4.*
   mytgl[24] + mytgl[39] + 5./4.*mytgl[17] - 763./32. - 1./3.*mytgl[46]
    - 3./2.*mytgl[45] + 1./2.*mytgl[43] - 3*mytgl[42];
   mytgl[40]=1./4.*MMH;
   mytgl[39]=mytgl[39]*mytgl[40];
   mytgl[42]=pow(mytgl[13],2);
   mytgl[43]=mytgl[4]*mytgl[42];
   mytgl[43]=mytgl[43] - mytgl[5];
   mytgl[45]=mytgl[11]*mytgl[13];
   mytgl[47]=mytgl[13] - 1./4.*mytgl[26];
   mytgl[48]=3./8.*mytgl[15];
   mytgl[49]=mytgl[13]*mytgl[48];
   mytgl[39]=mytgl[39] - 1./2.*mytgl[6] + mytgl[49] + 37./8.*mytgl[20]
    + 29./16.*mytgl[19] + 3*mytgl[47] - 7./4.*mytgl[45] + 3./16.*
   mytgl[43];
   mytgl[39]=MMH*mytgl[39];
   mytgl[43]=S2*pow(MMH,2);
   mytgl[39]=243./16.*mytgl[43] + 25./4.*mytgl[42] + mytgl[39];
   mytgl[43]=3*mytgl[11];
   mytgl[47]= - 3 + mytgl[14];
   mytgl[47]=5./2.*mytgl[47] + mytgl[11];
   mytgl[47]=mytgl[47]*mytgl[43];
   mytgl[49]=1./8.*mytgl[16];
   mytgl[50]=3./4.*mytgl[16];
   mytgl[51]=mytgl[50] - 39./4. + 11*mytgl[11];
   mytgl[51]=mytgl[51]*mytgl[49];
   mytgl[52]=5./8.*mytgl[34];
   mytgl[53]=3*mytgl[14];
   mytgl[47]=1./4.*mytgl[37] + 13./8.*mytgl[18] - mytgl[24] + 75./16.*
   mytgl[28] + mytgl[51] + mytgl[47] - mytgl[53] - 3./32.*mytgl[17] - 
   mytgl[52] - 657./128. - mytgl[46];
   mytgl[51]= - 1 - mytgl[14];
   mytgl[43]=mytgl[51]*mytgl[43];
   mytgl[51]= - 111 - 29./2.*mytgl[46];
   mytgl[43]=mytgl[43] + 1./16.*mytgl[51] + mytgl[53];
   mytgl[43]=mytgl[2]*mytgl[43];
   mytgl[43]=mytgl[23] + mytgl[43];
   mytgl[43]=MMt*mytgl[43];
   mytgl[43]=mytgl[43] + mytgl[25];
   mytgl[51]=7./8. - mytgl[14];
   mytgl[51]=mytgl[13]*mytgl[51];
   mytgl[51]= - mytgl[45] + mytgl[51];
   mytgl[51]=3*mytgl[51] - 169./16.*mytgl[19];
   mytgl[51]= - 9./8.*mytgl[6] + 3./2.*mytgl[9] + 1./2.*mytgl[51] - 6*
   mytgl[20];
   mytgl[51]=mytgl[2]*mytgl[51];
   mytgl[54]=mytgl[32] + mytgl[30];
   mytgl[54]=3*mytgl[22] + 1./4.*mytgl[54] - mytgl[23];
   mytgl[40]=mytgl[54]*mytgl[40];
   mytgl[40]=mytgl[40] - 115./32.*mytgl[21] + 1./4.*mytgl[47] + 
   mytgl[51] + 1./2.*mytgl[43];
   mytgl[40]=MMt*mytgl[40];
   mytgl[43]=3./4.*mytgl[14] - 3./16.*mytgl[17] - 9./4.*mytgl[27] + 
   mytgl[52] + 17 + 9./8.*mytgl[46];
   mytgl[46]=mytgl[11] - 1;
   mytgl[47]=mytgl[46]*mytgl[48];
   mytgl[48]=1./8.*MMH;
   mytgl[51]=mytgl[33] - 9*mytgl[22] + mytgl[31] + mytgl[29] - 
   mytgl[23];
   mytgl[51]=mytgl[51]*mytgl[48];
   mytgl[52]=39./4. - mytgl[14];
   mytgl[52]=3./2.*mytgl[52] - mytgl[11];
   mytgl[52]=mytgl[11]*mytgl[52];
   mytgl[43]=mytgl[51] - 5./4.*mytgl[25] + 35./4.*mytgl[21] + 3./16.*
   mytgl[37] + 35./16.*mytgl[18] + mytgl[47] + 5./8.*mytgl[28] - 11./8.
   *mytgl[35] + 1./32.*mytgl[16] + 1./2.*mytgl[43] + mytgl[52];
   mytgl[43]=mytgl[43]*mytgl[41];
   mytgl[47]= - mytgl[50] - 431./16. + mytgl[53];
   mytgl[47]=mytgl[13]*mytgl[47];
   mytgl[45]=129./8.*mytgl[19] + 9*mytgl[45] - 15./8.*mytgl[26] + 
   mytgl[47];
   mytgl[42]=mytgl[2]*mytgl[42];
   mytgl[42]=mytgl[43] + 15./8.*mytgl[6] - 93./8.*mytgl[42] + 5./16.*
   mytgl[9] + 1./4.*mytgl[45] + mytgl[20];
   mytgl[40]=1./2.*mytgl[42] + mytgl[40];
   mytgl[40]=MMt*mytgl[40];
   mytgl[42]=MMH*mytgl[4];
   mytgl[43]=203./2. + mytgl[11];
   mytgl[45]= - mytgl[4]*mytgl[13];
   mytgl[43]=mytgl[42] + mytgl[45] + 1./2.*mytgl[43] - 3*mytgl[15];
   mytgl[43]=mytgl[43]*mytgl[48];
   mytgl[43]= - 3*mytgl[13] + mytgl[43];
   mytgl[38]= - mytgl[38] - 23./16. - mytgl[14];
   mytgl[45]=mytgl[2]*mytgl[13];
   mytgl[38]=3*mytgl[38] + 109./4.*mytgl[45];
   mytgl[45]= - 1./8. + mytgl[14];
   mytgl[45]=mytgl[49] + 3*mytgl[45] + mytgl[11];
   mytgl[47]=MMt*mytgl[2];
   mytgl[45]=mytgl[45]*mytgl[47];
   mytgl[38]=1./4.*mytgl[38] + mytgl[45];
   mytgl[38]=MMt*mytgl[38];
   mytgl[45]=mytgl[10]*MMH;
   mytgl[48]=MMt*pow(mytgl[2],2);
   mytgl[48]=5./2.*mytgl[2] - 7*mytgl[48];
   mytgl[48]=MMt*mytgl[48];
   mytgl[48]= - 25./2. + 9*mytgl[48];
   mytgl[48]=mytgl[12]*mytgl[48];
   mytgl[38]=1./16.*mytgl[48] - 9./16.*mytgl[45] + 1./2.*mytgl[43] + 
   mytgl[38];
   mytgl[38]=mytgl[12]*mytgl[38];
   mytgl[43]= - MMH*mytgl[44];
   mytgl[44]=MMt*mytgl[46];
   mytgl[43]=mytgl[44] + mytgl[13] + mytgl[43];
   mytgl[43]=mytgl[43]*mytgl[45];
   mytgl[44]=3./2. - mytgl[47];
   mytgl[44]=MMt*mytgl[44];
   mytgl[41]= - mytgl[41] + mytgl[44];
   mytgl[41]=mytgl[7]*mytgl[41];
   mytgl[42]= - 5./2. - mytgl[42];
   mytgl[42]=mytgl[8]*MMH*mytgl[42];
   mytgl[38]=1./2.*mytgl[38] + 9./16.*mytgl[41] + 9./32.*mytgl[43] + 1./
   32.*mytgl[42] + 1./4.*mytgl[39] + mytgl[40];

      return mytgl[38]*pow(mytgl[3],2)*pow(mytgl[1],4);
}
