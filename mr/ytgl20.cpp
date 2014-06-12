#include <tt.hpp>
std::complex<long double> tt::mygl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[65];

    mytgl[1]=pow(SW,-1);
    mytgl[2]=pow(MMH,-1);
    mytgl[3]=pow(MMW,-1);
    mytgl[4]=pow(MMt,-1);
    mytgl[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    mytgl[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    mytgl[7]=Tsil::B(MMH,MMt,MMt,mu2);
    mytgl[8]=Tsil::B(MMt,MMt,MMH,mu2);
    mytgl[9]=log(pow(mu2,-1)*MMt);
    mytgl[10]=log(pow(mu2,-1)*MMH);
    mytgl[11]=Tsil::B(0,MMH,MMt,mu2);
    mytgl[12]=Tsil::B(0,0,MMH,mu2);
    mytgl[13]=Tsil::B(0,0,MMt,mu2);
    mytgl[14]=Tsil::Beps(MMH,MMt,MMt,mu2);
    mytgl[15]=log(MMH);
    mytgl[16]=log(MMt);
    mytgl[17]=prot0ttHt->Tuxv(0);
    mytgl[18]=protHHttH->M(0);
    mytgl[19]=protHttHt->M(0);
    mytgl[20]=protHHttH->Uxzuv(0);
    mytgl[21]=protHttHt->Txuv(0);
    mytgl[22]=protHHttH->Suxv(0);
    mytgl[23]=protHHttH->Uuyxv(0);
    mytgl[24]=protWt000->Tyzv(0);
    mytgl[25]=protH0tt0->M(0);
    mytgl[26]=protH0t00->M(0);
    mytgl[27]=prot0Htt0->M(0);
    mytgl[28]=prot0H0t0->M(0);
    mytgl[29]=prot00ttH->M(0);
    mytgl[30]=protH0tt0->Uzxyv(0);
    mytgl[31]=protH0t00->Uxzuv(0);
    mytgl[32]=protH0t00->Uuyxv(0);
    mytgl[33]=protH0t00->Txuv(0);
    mytgl[34]=1/(MMt - MMH);
    mytgl[35]=Tsil::I2(0,MMH,MMt,mu2);
    mytgl[36]=1/(4*MMt - MMH);
   mytgl[37]=1./2.*mytgl[15];
   mytgl[38]=mytgl[37]*mytgl[36];
   mytgl[39]=mytgl[16] - 3;
   mytgl[40]=mytgl[39]*mytgl[36];
   mytgl[38]=mytgl[38] + mytgl[40];
   mytgl[41]=mytgl[38]*mytgl[10];
   mytgl[42]=mytgl[9]*mytgl[36];
   mytgl[43]=1./2.*mytgl[16];
   mytgl[44]=1 - mytgl[43];
   mytgl[44]=mytgl[44]*mytgl[42];
   mytgl[44]=mytgl[41] + 3*mytgl[44] - 1./2.*mytgl[18];
   mytgl[45]=mytgl[27] + mytgl[25];
   mytgl[46]=1./2.*mytgl[19];
   mytgl[44]= - mytgl[46] + 9*mytgl[44] + 1./2.*mytgl[45];
   mytgl[45]=pow(Pi,2);
   mytgl[47]=1./16.*mytgl[45];
   mytgl[48]=7./4.*mytgl[16];
   mytgl[49]= - 3 + mytgl[48];
   mytgl[49]=mytgl[9]*mytgl[49];
   mytgl[49]= - 21./2. + mytgl[49];
   mytgl[50]=1 - 1./8.*mytgl[16];
   mytgl[50]=3*mytgl[50] - 13./16.*mytgl[15];
   mytgl[50]=mytgl[10]*mytgl[50];
   mytgl[49]= - mytgl[47] + 1./4.*mytgl[49] + mytgl[50];
   mytgl[49]=mytgl[34]*mytgl[49];
   mytgl[44]=1./4.*mytgl[29] + 1./2.*mytgl[44] + 9*mytgl[49];
   mytgl[49]=1./2.*MMH;
   mytgl[44]=mytgl[49]*mytgl[44];
   mytgl[50]=1./4.*mytgl[9];
   mytgl[51]=35./4. - mytgl[7];
   mytgl[51]=7*mytgl[51] - 25./4.*mytgl[16];
   mytgl[51]=mytgl[51]*mytgl[50];
   mytgl[52]=3./2.*mytgl[10];
   mytgl[53]=5*mytgl[7];
   mytgl[54]=239./4. - mytgl[53];
   mytgl[54]= - 55./16.*mytgl[15] + 1./4.*mytgl[54] - 5*mytgl[16];
   mytgl[54]=mytgl[54]*mytgl[52];
   mytgl[55]=1./2.*mytgl[9];
   mytgl[56]=mytgl[55] - 1 + 1./2.*mytgl[7];
   mytgl[57]=3./4.*mytgl[12];
   mytgl[58]=mytgl[56]*mytgl[57];
   mytgl[59]=mytgl[34]*mytgl[35];
   mytgl[60]=mytgl[36]*mytgl[6];
   mytgl[61]= - 205./4. + mytgl[24];
   mytgl[62]=301./16. - mytgl[7];
   mytgl[62]=mytgl[7]*mytgl[62];
   mytgl[63]=mytgl[7] + mytgl[10] - 5./4.;
   mytgl[64]=mytgl[8]*mytgl[63];
   mytgl[44]=mytgl[44] + 35./4.*mytgl[17] + mytgl[58] - 3./32.*
   mytgl[11] - 9./2.*mytgl[59] - 3./2.*mytgl[64] + 1./32.*mytgl[13] - 
   203./192.*mytgl[45] + mytgl[54] + mytgl[51] - 5./4.*mytgl[21] + 35./
   16.*mytgl[14] - 9./8.*mytgl[23] - 11./8.*mytgl[31] + 3./16.*
   mytgl[33] + 5./16.*mytgl[30] + 9./4.*mytgl[60] + 5./8.*mytgl[61] + 
   mytgl[62];
   mytgl[51]=pow(mytgl[3],2);
   mytgl[44]=MMH*mytgl[51]*mytgl[44];
   mytgl[54]=mytgl[6] - 5./8.*mytgl[22];
   mytgl[54]=mytgl[54]*mytgl[51];
   mytgl[44]=3./2.*mytgl[54] + mytgl[44];
   mytgl[54]=1679./8. + 15*mytgl[24];
   mytgl[58]= - 29 + 3*mytgl[7];
   mytgl[58]=mytgl[7]*mytgl[58];
   mytgl[61]=mytgl[2]*mytgl[6];
   mytgl[54]=mytgl[45] + 9*mytgl[61] + 13./8.*mytgl[14] + 1./4.*
   mytgl[33] - 5./8.*mytgl[30] + 45./2.*mytgl[60] + mytgl[58] + 5./16.*
   mytgl[54] - mytgl[20];
   mytgl[58]= - 1 - mytgl[43];
   mytgl[58]=mytgl[58]*mytgl[42];
   mytgl[61]=mytgl[28] + mytgl[26];
   mytgl[58]=27*mytgl[58] + 1./2.*mytgl[61] + 189*mytgl[36];
   mytgl[61]=3 - 17./16.*mytgl[16];
   mytgl[61]=mytgl[9]*mytgl[61];
   mytgl[62]=19./4.*mytgl[15] - 15 + mytgl[43];
   mytgl[62]=mytgl[10]*mytgl[62];
   mytgl[61]=1./4.*mytgl[62] + 7./8. + mytgl[61];
   mytgl[47]=3*mytgl[61] + mytgl[47];
   mytgl[61]=3./2.*mytgl[34];
   mytgl[47]=mytgl[47]*mytgl[61];
   mytgl[62]=mytgl[45]*mytgl[36];
   mytgl[41]=mytgl[47] + 9./4.*mytgl[62] - mytgl[19] + 45./2.*mytgl[41]
    + 1./2.*mytgl[58] + 3*mytgl[18];
   mytgl[41]=mytgl[41]*mytgl[49];
   mytgl[47]= - 43./16. + mytgl[7];
   mytgl[47]=9*mytgl[47] - 311./16.*mytgl[16];
   mytgl[47]=mytgl[47]*mytgl[50];
   mytgl[50]=3*mytgl[8];
   mytgl[53]=mytgl[9] - 7 + mytgl[53];
   mytgl[53]=1./4.*mytgl[53] + mytgl[10];
   mytgl[53]=mytgl[53]*mytgl[50];
   mytgl[58]=383./64.*mytgl[15] + 425./32.*mytgl[16] - 67./2. + 
   mytgl[7];
   mytgl[58]=mytgl[10]*mytgl[58];
   mytgl[62]=3*mytgl[9];
   mytgl[64]= - mytgl[62] - 35./4. + 11*mytgl[7];
   mytgl[64]=3./8.*mytgl[13] + 1./2.*mytgl[64] + mytgl[10];
   mytgl[64]=mytgl[13]*mytgl[64];
   mytgl[41]=mytgl[41] - 115./16.*mytgl[17] - 3./64.*mytgl[11] + 45./16.
   *mytgl[59] + mytgl[53] + 1./8.*mytgl[64] + mytgl[58] + mytgl[47] + 
   mytgl[21] + 1./2.*mytgl[54];
   mytgl[39]= - mytgl[2]*mytgl[39]*mytgl[42];
   mytgl[47]=mytgl[2]*mytgl[36];
   mytgl[39]= - 7./2.*mytgl[47] + mytgl[39];
   mytgl[53]=1./4.*mytgl[45];
   mytgl[47]= - mytgl[53]*mytgl[47];
   mytgl[39]=3*mytgl[39] + mytgl[47];
   mytgl[47]=1./4.*mytgl[16];
   mytgl[54]=1 - mytgl[47];
   mytgl[54]=mytgl[54]*mytgl[62];
   mytgl[58]=mytgl[37] - mytgl[16];
   mytgl[62]=mytgl[10]*mytgl[58];
   mytgl[54]=1./2.*mytgl[62] - 7./2. + mytgl[54];
   mytgl[54]=3*mytgl[54] - mytgl[53];
   mytgl[54]=mytgl[34]*mytgl[2]*mytgl[54];
   mytgl[39]=3*mytgl[39] + 1./16.*mytgl[54];
   mytgl[39]=MMt*mytgl[39];
   mytgl[54]=3*mytgl[36];
   mytgl[62]= - 15./4. + 2*mytgl[16];
   mytgl[62]=mytgl[62]*mytgl[54];
   mytgl[64]=61./32.*mytgl[16] - 7 - mytgl[7];
   mytgl[64]=mytgl[2]*mytgl[64];
   mytgl[62]=mytgl[62] + 1./2.*mytgl[64];
   mytgl[62]=mytgl[9]*mytgl[62];
   mytgl[58]=mytgl[2]*mytgl[58];
   mytgl[38]= - 3*mytgl[38] - 1./16.*mytgl[58];
   mytgl[38]=mytgl[38]*mytgl[52];
   mytgl[52]=7./2. - 3*mytgl[60];
   mytgl[52]=mytgl[2]*mytgl[52];
   mytgl[52]= - 7./4.*mytgl[36] + mytgl[52];
   mytgl[38]=mytgl[39] + mytgl[38] + 3./2.*mytgl[52] + mytgl[62];
   mytgl[39]=1./8.*mytgl[45];
   mytgl[52]= - 15 + 17./4.*mytgl[16];
   mytgl[52]=mytgl[9]*mytgl[52];
   mytgl[52]=21./2. + mytgl[52];
   mytgl[47]=1 + mytgl[47];
   mytgl[47]=3*mytgl[47] - 11./8.*mytgl[15];
   mytgl[47]=mytgl[10]*mytgl[47];
   mytgl[58]= - mytgl[35]*mytgl[2];
   mytgl[47]=mytgl[39] + mytgl[58] + 1./2.*mytgl[52] + mytgl[47];
   mytgl[47]=mytgl[34]*mytgl[47];
   mytgl[52]= - mytgl[36] + 3./2.*mytgl[2];
   mytgl[52]=mytgl[52]*mytgl[45];
   mytgl[50]= - mytgl[50]*mytgl[2]*mytgl[56];
   mytgl[38]=9./16.*mytgl[47] + mytgl[50] + 3./16.*mytgl[52] + 
   mytgl[46] + 3*mytgl[38];
   mytgl[38]=MMt*mytgl[38];
   mytgl[38]=1./2.*mytgl[41] + mytgl[38];
   mytgl[38]=MMt*mytgl[51]*mytgl[38];
   mytgl[38]=1./4.*mytgl[44] + mytgl[38];
   mytgl[38]=MMt*mytgl[38];
   mytgl[41]= - 21*mytgl[36] + mytgl[4];
   mytgl[44]= - mytgl[33]*mytgl[4];
   mytgl[42]=mytgl[16]*mytgl[42];
   mytgl[46]=mytgl[11]*mytgl[4];
   mytgl[41]= - mytgl[46] + 9./2.*mytgl[42] + 3./2.*mytgl[41] + 
   mytgl[44];
   mytgl[42]= - 9./2.*mytgl[36] + mytgl[4];
   mytgl[37]=mytgl[42]*mytgl[37];
   mytgl[37]=mytgl[37] - 9./2.*mytgl[40] - mytgl[4];
   mytgl[37]=mytgl[10]*mytgl[37];
   mytgl[40]= - mytgl[54] + mytgl[4];
   mytgl[39]=mytgl[40]*mytgl[39];
   mytgl[40]= - mytgl[9]*mytgl[43];
   mytgl[40]=7 + mytgl[40];
   mytgl[42]=3./4.*mytgl[15] - 3 + mytgl[43];
   mytgl[42]=mytgl[10]*mytgl[42];
   mytgl[40]=1./2.*mytgl[40] + mytgl[42];
   mytgl[40]=3*mytgl[40] + mytgl[53];
   mytgl[40]=mytgl[40]*mytgl[61];
   mytgl[37]=mytgl[40] + mytgl[39] + mytgl[37] + 1./2.*mytgl[41];
   mytgl[37]=mytgl[37]*mytgl[49];
   mytgl[39]=3*mytgl[20] - 1645./8. + 243*S2;
   mytgl[40]= - 117./4. + mytgl[7];
   mytgl[40]=mytgl[7]*mytgl[40];
   mytgl[39]=3*mytgl[31] - 9*mytgl[60] + 1./2.*mytgl[39] + mytgl[40];
   mytgl[40]=3 + mytgl[48];
   mytgl[40]=mytgl[40]*mytgl[55];
   mytgl[41]=mytgl[21] - mytgl[23];
   mytgl[42]= - 11./2.*mytgl[16] + 59 + 19./2.*mytgl[7];
   mytgl[42]=1./2.*mytgl[42] - 3*mytgl[15];
   mytgl[42]=mytgl[10]*mytgl[42];
   mytgl[39]=mytgl[42] + mytgl[40] - 21./4.*mytgl[14] + 1./2.*mytgl[39]
    - 3*mytgl[41];
   mytgl[40]= - mytgl[63]*mytgl[57];
   mytgl[37]=mytgl[37] - 23./8.*mytgl[17] + mytgl[40] + 1./4.*mytgl[32]
    + 5./8.*mytgl[11] + 9./4.*mytgl[59] + 1./2.*mytgl[39] - 1./3.*
   mytgl[45];
   mytgl[37]=MMH*mytgl[37];
   mytgl[39]= - 3./4.*mytgl[5] - 13./2.*mytgl[6] - 3*mytgl[22];
   mytgl[37]=1./2.*mytgl[39] + mytgl[37];
   mytgl[37]=MMH*mytgl[51]*mytgl[37];
   mytgl[37]=1./8.*mytgl[37] + mytgl[38];

      return mytgl[37]*pow(mytgl[1],4);
}
