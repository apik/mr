#include <tt.hpp>
std::complex<long double> tt::mytgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mytgl[68];

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
    mytgl[17]=Tfin(MMH,MMt,MMH);
    mytgl[18]=Tfin(MMH,MMt,0);
    mytgl[19]=Tfin(MMH,0,0);
    mytgl[20]=Tfin(MMt,0,0);
    mytgl[21]=Mfin(MMH,MMH,MMt,MMt,MMH);
    mytgl[22]=Mfin(MMH,MMt,MMt,MMH,MMt);
    mytgl[23]=Mfin(MMH,0,MMt,MMt,0);
    mytgl[24]=Mfin(MMH,0,MMt,0,0);
    mytgl[25]=Mfin(0,MMH,MMt,MMt,0);
    mytgl[26]=Mfin(0,MMH,0,MMt,0);
    mytgl[27]=Mfin(0,0,MMt,MMt,MMH);
    mytgl[28]=Sfin(MMH,MMH,MMt);
    mytgl[29]=Ufin(MMH,MMt,MMH,MMt);
    mytgl[30]=Ufin(MMH,MMt,0,0);
    mytgl[31]=Ufin(MMt,MMH,MMH,MMH);
    mytgl[32]=Ufin(MMt,MMH,0,0);
    mytgl[33]=Ufin(MMt,0,0,MMH);
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
   mytgl[45]=mytgl[25] + mytgl[27] - 9*mytgl[21] - mytgl[22] + 
   mytgl[23];
   mytgl[44]=9*mytgl[41] + 27*mytgl[44] + 1./2.*mytgl[45];
   mytgl[45]=pow(Pi,2);
   mytgl[46]=1./16.*mytgl[45];
   mytgl[47]=7./4.*mytgl[16];
   mytgl[48]= - 3 + mytgl[47];
   mytgl[48]=mytgl[9]*mytgl[48];
   mytgl[48]= - 21./2. + mytgl[48];
   mytgl[49]=1 - 1./8.*mytgl[16];
   mytgl[49]=3*mytgl[49] - 13./16.*mytgl[15];
   mytgl[49]=mytgl[10]*mytgl[49];
   mytgl[48]= - mytgl[46] + 1./4.*mytgl[48] + mytgl[49];
   mytgl[48]=mytgl[34]*mytgl[48];
   mytgl[44]=1./2.*mytgl[44] + 9*mytgl[48];
   mytgl[48]=1./2.*MMH;
   mytgl[44]=mytgl[44]*mytgl[48];
   mytgl[49]=5*mytgl[7];
   mytgl[50]=239./4. - mytgl[49];
   mytgl[50]= - 55./16.*mytgl[15] + 1./4.*mytgl[50] - 5*mytgl[16];
   mytgl[50]=mytgl[10]*mytgl[50];
   mytgl[51]=mytgl[7] + mytgl[10] - 5./4.;
   mytgl[52]=mytgl[8]*mytgl[51];
   mytgl[50]=mytgl[50] - mytgl[52];
   mytgl[52]=1./2.*mytgl[9];
   mytgl[53]=mytgl[52] - 1 + 1./2.*mytgl[7];
   mytgl[54]=3./4.*mytgl[12];
   mytgl[55]=mytgl[53]*mytgl[54];
   mytgl[56]=mytgl[34]*mytgl[35];
   mytgl[57]=mytgl[36]*mytgl[6];
   mytgl[58]= - 205./4. + mytgl[20];
   mytgl[59]=301./16. - mytgl[7];
   mytgl[59]=mytgl[7]*mytgl[59];
   mytgl[60]=35./4. - mytgl[7];
   mytgl[60]=7*mytgl[60] - 25./4.*mytgl[16];
   mytgl[60]=mytgl[9]*mytgl[60];
   mytgl[44]=mytgl[44] - 5./4.*mytgl[17] + mytgl[55] - 11./8.*mytgl[32]
    - 3./32.*mytgl[11] - 9./2.*mytgl[56] + 1./32.*mytgl[13] - 203./192.
   *mytgl[45] + 3./16.*mytgl[19] + 35./4.*mytgl[18] + 1./4.*mytgl[60]
    + 35./16.*mytgl[14] - 9./8.*mytgl[31] + 5./16.*mytgl[30] + 9./4.*
   mytgl[57] + 5./8.*mytgl[58] + mytgl[59] + 3./2.*mytgl[50];
   mytgl[50]=pow(mytgl[3],2);
   mytgl[55]=mytgl[50]*MMH;
   mytgl[44]=mytgl[44]*mytgl[55];
   mytgl[58]= - 5./8.*mytgl[28] + mytgl[6];
   mytgl[58]=mytgl[58]*mytgl[50];
   mytgl[44]=3./2.*mytgl[58] + mytgl[44];
   mytgl[58]=1./8.*mytgl[45];
   mytgl[59]= - 15 + 17./4.*mytgl[16];
   mytgl[59]=mytgl[9]*mytgl[59];
   mytgl[59]=21./2. + mytgl[59];
   mytgl[60]=1./4.*mytgl[16];
   mytgl[61]=1 + mytgl[60];
   mytgl[61]=3*mytgl[61] - 11./8.*mytgl[15];
   mytgl[61]=mytgl[10]*mytgl[61];
   mytgl[62]= - mytgl[35]*mytgl[2];
   mytgl[59]=mytgl[58] + mytgl[62] + 1./2.*mytgl[59] + mytgl[61];
   mytgl[59]=mytgl[34]*mytgl[59];
   mytgl[39]= - mytgl[2]*mytgl[39]*mytgl[42];
   mytgl[61]=mytgl[2]*mytgl[36];
   mytgl[39]= - 7./2.*mytgl[61] + mytgl[39];
   mytgl[62]=1./4.*mytgl[45];
   mytgl[61]= - mytgl[62]*mytgl[61];
   mytgl[39]=3*mytgl[39] + mytgl[61];
   mytgl[61]=3*mytgl[9];
   mytgl[60]=1 - mytgl[60];
   mytgl[60]=mytgl[60]*mytgl[61];
   mytgl[63]=mytgl[37] - mytgl[16];
   mytgl[64]=mytgl[10]*mytgl[63];
   mytgl[60]=1./2.*mytgl[64] - 7./2. + mytgl[60];
   mytgl[60]=3*mytgl[60] - mytgl[62];
   mytgl[60]=mytgl[34]*mytgl[2]*mytgl[60];
   mytgl[39]=3*mytgl[39] + 1./16.*mytgl[60];
   mytgl[39]=MMt*mytgl[39];
   mytgl[60]=9*mytgl[2];
   mytgl[64]=7./2. - 3*mytgl[57];
   mytgl[64]=mytgl[64]*mytgl[60];
   mytgl[64]=mytgl[64] - 63./4.*mytgl[36] + mytgl[22];
   mytgl[65]=3*mytgl[36];
   mytgl[66]= - 15./4. + 2*mytgl[16];
   mytgl[66]=mytgl[66]*mytgl[65];
   mytgl[67]=61./32.*mytgl[16] - 7 - mytgl[7];
   mytgl[67]=mytgl[2]*mytgl[67];
   mytgl[66]=mytgl[66] + 1./2.*mytgl[67];
   mytgl[66]=mytgl[66]*mytgl[61];
   mytgl[63]=mytgl[2]*mytgl[63];
   mytgl[38]= - 3*mytgl[38] - 1./16.*mytgl[63];
   mytgl[38]=mytgl[10]*mytgl[38];
   mytgl[63]= - mytgl[36] + 3./2.*mytgl[2];
   mytgl[63]=mytgl[63]*mytgl[45];
   mytgl[67]=3*mytgl[8];
   mytgl[53]= - mytgl[67]*mytgl[2]*mytgl[53];
   mytgl[38]=3*mytgl[39] + 9./16.*mytgl[59] + mytgl[53] + 3./16.*
   mytgl[63] + 9./2.*mytgl[38] + 1./2.*mytgl[64] + mytgl[66];
   mytgl[38]=MMt*mytgl[38];
   mytgl[39]= - 43./16. + mytgl[7];
   mytgl[39]=9*mytgl[39] - 311./16.*mytgl[16];
   mytgl[39]=mytgl[39]*mytgl[52];
   mytgl[53]=1679./8. + 15*mytgl[20];
   mytgl[59]= - 29 + 3*mytgl[7];
   mytgl[59]=mytgl[7]*mytgl[59];
   mytgl[60]=mytgl[6]*mytgl[60];
   mytgl[39]= - mytgl[29] + mytgl[45] - 115./8.*mytgl[18] + mytgl[39]
    + mytgl[60] + 13./8.*mytgl[14] - 5./8.*mytgl[30] + 45./2.*mytgl[57]
    + 5./16.*mytgl[53] + mytgl[59];
   mytgl[53]= - 1 - mytgl[43];
   mytgl[53]=mytgl[53]*mytgl[42];
   mytgl[59]=3 - 17./16.*mytgl[16];
   mytgl[59]=mytgl[9]*mytgl[59];
   mytgl[60]=19./4.*mytgl[15] - 15 + mytgl[43];
   mytgl[60]=mytgl[10]*mytgl[60];
   mytgl[59]=1./4.*mytgl[60] + 7./8. + mytgl[59];
   mytgl[46]=3*mytgl[59] + mytgl[46];
   mytgl[59]=3./2.*mytgl[34];
   mytgl[46]=mytgl[46]*mytgl[59];
   mytgl[60]=mytgl[45]*mytgl[36];
   mytgl[63]=mytgl[24] + mytgl[26];
   mytgl[63]=1./2.*mytgl[63] + 189*mytgl[36];
   mytgl[41]=mytgl[46] + 9./4.*mytgl[60] + 45./2.*mytgl[41] + 27./2.*
   mytgl[53] + 3*mytgl[21] + 1./2.*mytgl[63] - mytgl[22];
   mytgl[41]=mytgl[41]*mytgl[48];
   mytgl[46]=mytgl[9] - 7 + mytgl[49];
   mytgl[46]=1./4.*mytgl[46] + mytgl[10];
   mytgl[46]=mytgl[46]*mytgl[67];
   mytgl[49]= - mytgl[61] - 35./4. + 11*mytgl[7];
   mytgl[49]=3./8.*mytgl[13] + 1./2.*mytgl[49] + mytgl[10];
   mytgl[49]=mytgl[13]*mytgl[49];
   mytgl[49]=mytgl[49] + mytgl[19];
   mytgl[53]=383./64.*mytgl[15] + 425./32.*mytgl[16] - 67./2. + 
   mytgl[7];
   mytgl[53]=mytgl[10]*mytgl[53];
   mytgl[39]=mytgl[41] + mytgl[17] - 3./64.*mytgl[11] + 45./16.*
   mytgl[56] + mytgl[46] + mytgl[53] + 1./8.*mytgl[49] + 1./2.*
   mytgl[39];
   mytgl[38]=1./2.*mytgl[39] + mytgl[38];
   mytgl[38]=MMt*mytgl[50]*mytgl[38];
   mytgl[38]=1./4.*mytgl[44] + mytgl[38];
   mytgl[38]=MMt*mytgl[38];
   mytgl[39]= - 9./2.*mytgl[36] + mytgl[4];
   mytgl[37]=mytgl[39]*mytgl[37];
   mytgl[37]=mytgl[37] - 9./2.*mytgl[40] - mytgl[4];
   mytgl[37]=mytgl[10]*mytgl[37];
   mytgl[39]= - mytgl[65] + mytgl[4];
   mytgl[39]=mytgl[39]*mytgl[58];
   mytgl[40]= - mytgl[16]*mytgl[52];
   mytgl[40]=7 + mytgl[40];
   mytgl[41]=3./4.*mytgl[15] - 3 + mytgl[43];
   mytgl[41]=mytgl[10]*mytgl[41];
   mytgl[40]=1./2.*mytgl[40] + mytgl[41];
   mytgl[40]=3*mytgl[40] + mytgl[62];
   mytgl[40]=mytgl[40]*mytgl[59];
   mytgl[41]=mytgl[11] + mytgl[19];
   mytgl[41]= - 1./2.*mytgl[41];
   mytgl[41]=mytgl[4]*mytgl[41];
   mytgl[42]=mytgl[16]*mytgl[42];
   mytgl[42]=3*mytgl[42] - 21*mytgl[36] + mytgl[4];
   mytgl[37]=mytgl[40] + mytgl[39] + 3./4.*mytgl[42] + mytgl[37] + 
   mytgl[41];
   mytgl[37]=mytgl[37]*mytgl[48];
   mytgl[39]= - 1645./8. + 243*S2;
   mytgl[40]= - 117./4. + mytgl[7];
   mytgl[40]=mytgl[7]*mytgl[40];
   mytgl[39]=mytgl[33] - 9*mytgl[57] + 1./2.*mytgl[39] + mytgl[40];
   mytgl[40]=3 + mytgl[47];
   mytgl[40]=mytgl[40]*mytgl[52];
   mytgl[41]= - 11./2.*mytgl[16] + 59 + 19./2.*mytgl[7];
   mytgl[41]=1./2.*mytgl[41] - 3*mytgl[15];
   mytgl[41]=mytgl[10]*mytgl[41];
   mytgl[39]=mytgl[41] - 23./4.*mytgl[18] + mytgl[40] - 21./4.*
   mytgl[14] + 1./2.*mytgl[39] + 3*mytgl[31];
   mytgl[40]= - mytgl[51]*mytgl[54];
   mytgl[37]=mytgl[37] + 3./8.*mytgl[29] - 3./2.*mytgl[17] + mytgl[40]
    + 3./4.*mytgl[32] + 5./8.*mytgl[11] + 9./4.*mytgl[56] + 1./2.*
   mytgl[39] - 1./3.*mytgl[45];
   mytgl[37]=mytgl[37]*mytgl[55];
   mytgl[39]= - 3./4.*mytgl[5] - 3*mytgl[28] - 13./2.*mytgl[6];
   mytgl[39]=mytgl[39]*mytgl[50];
   mytgl[37]=1./2.*mytgl[39] + mytgl[37];
   mytgl[37]=MMH*mytgl[37];
   mytgl[37]=1./8.*mytgl[37] + mytgl[38];

      return mytgl[37]*pow(mytgl[1],4);
}
