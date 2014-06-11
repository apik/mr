#include <HH.hpp>
std::complex<long double> HH::lamgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mlamgl[53];

    mlamgl[1]=pow(SW,-1);
    mlamgl[2]=pow(MMH,-1);
    mlamgl[3]=pow(MMW,-1);
    mlamgl[4]=Tsil::I2(MMH,MMH,MMH,mu2);
    mlamgl[5]=Tsil::I2(MMH,MMt,MMt,mu2);
    mlamgl[6]=Tsil::I2(0,MMH,MMt,mu2);
    mlamgl[7]=Tsil::B(MMH,MMH,MMH,mu2);
    mlamgl[8]=Tsil::B(MMt,MMt,MMH,mu2);
    mlamgl[9]=log(pow(mu2,-1)*MMt);
    mlamgl[10]=Tsil::B(MMH,MMt,MMt,mu2);
    mlamgl[11]=log(pow(mu2,-1)*MMH);
    mlamgl[12]=Tsil::B(0,MMt,MMH,mu2);
    mlamgl[13]=Tsil::B(0,0,MMt,mu2);
    mlamgl[14]=Tsil::Beps(MMH,MMH,MMH,mu2);
    mlamgl[15]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mlamgl[16]=log(MMH);
    mlamgl[17]=log(MMt);
    mlamgl[18]=prottttt0->Suxv(0);
    mlamgl[19]=protHHHHH->M(0);
    mlamgl[20]=protHtHtt->M(0);
    mlamgl[21]=protttttH->M(0);
    mlamgl[22]=protHHHHH->Uzxyv(0);
    mlamgl[23]=protHtHtt->Uzxyv(0);
    mlamgl[24]=protHtHtt->Uuyxv(0);
    mlamgl[25]=protHtHtt->Tyzv(0);
    mlamgl[26]=protHtHtt->Svyz(0);
    mlamgl[27]=Mfin(0,MMH,0,MMH,0);
    mlamgl[28]=Mfin(0,MMt,0,MMt,MMt);
    mlamgl[29]=Mfin(0,MMt,0,MMt,0);
    mlamgl[30]=Mfin(0,0,0,0,MMH);
    mlamgl[31]=Tfin(MMH,0,0);
    mlamgl[32]=Tfin(MMt,0,0);
    mlamgl[33]=Ufin(MMH,MMH,0,0);
    mlamgl[34]=Ufin(MMt,MMt,0,0);
    mlamgl[35]=1/(4*MMt - MMH);
   mlamgl[36]=1./4.*mlamgl[8];
   mlamgl[37]= - 1 - 3./2.*mlamgl[7];
   mlamgl[37]= - mlamgl[36] + 3*mlamgl[37] + 1./2.*mlamgl[10];
   mlamgl[38]=1./2.*mlamgl[9];
   mlamgl[37]=mlamgl[37]*mlamgl[38];
   mlamgl[39]= - mlamgl[24] + mlamgl[15] - mlamgl[12];
   mlamgl[40]=mlamgl[23] - mlamgl[14];
   mlamgl[41]=1./2.*mlamgl[32];
   mlamgl[42]=43./8. + 9*mlamgl[7];
   mlamgl[43]=3*mlamgl[7];
   mlamgl[44]= - 1./16.*mlamgl[8] - 1./4.*mlamgl[10] + 11./8. + 
   mlamgl[43];
   mlamgl[44]=mlamgl[8]*mlamgl[44];
   mlamgl[45]=5./2.*mlamgl[17] - 1 + 3./4.*mlamgl[16];
   mlamgl[45]=1./4.*mlamgl[45] - mlamgl[8];
   mlamgl[45]=mlamgl[11]*mlamgl[45];
   mlamgl[37]=mlamgl[41] + mlamgl[45] + mlamgl[37] + mlamgl[44] + 1./8.
   *mlamgl[10] + 1./2.*mlamgl[42] - 3*mlamgl[40] + 3./4.*mlamgl[39];
   mlamgl[37]=MMH*mlamgl[37];
   mlamgl[37]= - 3./4.*mlamgl[5] + mlamgl[37] - mlamgl[18] + 9./8.*
   mlamgl[6];
   mlamgl[39]=pow(mlamgl[1],4);
   mlamgl[37]=mlamgl[39]*mlamgl[37];
   mlamgl[40]=55./2. - 63*mlamgl[7];
   mlamgl[40]= - mlamgl[36] + 1./2.*mlamgl[40] + 7*mlamgl[10];
   mlamgl[40]=mlamgl[8]*mlamgl[40];
   mlamgl[42]=9*mlamgl[14];
   mlamgl[40]=mlamgl[40] - 3*mlamgl[10] + 9./4.*mlamgl[12] + 9*
   mlamgl[23] - 23 - mlamgl[42];
   mlamgl[44]=1./2.*mlamgl[17];
   mlamgl[45]=mlamgl[44] + 37 - 27*mlamgl[7];
   mlamgl[45]= - 3./2.*mlamgl[8] + 1./2.*mlamgl[45] - mlamgl[10];
   mlamgl[45]=mlamgl[45]*mlamgl[38];
   mlamgl[46]=1./2.*mlamgl[11];
   mlamgl[47]= - 15./2.*mlamgl[8] - 9./2.*mlamgl[17] + 5./2. - 
   mlamgl[16];
   mlamgl[47]=mlamgl[47]*mlamgl[46];
   mlamgl[48]=1./2.*MMH;
   mlamgl[49]=1./2.*mlamgl[21] - mlamgl[28] + 9*mlamgl[20];
   mlamgl[48]=mlamgl[49]*mlamgl[48];
   mlamgl[49]=5./4.*mlamgl[34];
   mlamgl[50]= - mlamgl[9] - 1./2. + mlamgl[8];
   mlamgl[50]=mlamgl[13]*mlamgl[50];
   mlamgl[40]= - mlamgl[49] + 1./4.*mlamgl[50] + mlamgl[48] + 
   mlamgl[47] + 5./4.*mlamgl[15] + mlamgl[45] + 1./2.*mlamgl[40] - 3*
   mlamgl[25];
   mlamgl[40]=mlamgl[39]*mlamgl[40];
   mlamgl[45]=mlamgl[39]*mlamgl[2];
   mlamgl[47]= - 5*mlamgl[26] - mlamgl[18] - 3*mlamgl[6];
   mlamgl[47]=1./4.*mlamgl[47] + mlamgl[5];
   mlamgl[47]=mlamgl[47]*mlamgl[45];
   mlamgl[40]=mlamgl[47] + mlamgl[40];
   mlamgl[47]=mlamgl[10] + mlamgl[24];
   mlamgl[48]= - 5*mlamgl[10] - 7./8.*mlamgl[8];
   mlamgl[48]=mlamgl[8]*mlamgl[48];
   mlamgl[50]=13./4.*mlamgl[8] - mlamgl[10] + 29./16. + mlamgl[17];
   mlamgl[50]=mlamgl[9]*mlamgl[50];
   mlamgl[51]= - 2*mlamgl[8] - 4 + mlamgl[16];
   mlamgl[51]=mlamgl[11]*mlamgl[51];
   mlamgl[52]= - mlamgl[38] + 1 - 5./2.*mlamgl[8];
   mlamgl[52]=mlamgl[13]*mlamgl[52];
   mlamgl[47]=mlamgl[49] + 1./2.*mlamgl[52] - 3./8.*mlamgl[32] + 
   mlamgl[51] - 13./4.*mlamgl[15] + mlamgl[50] - 5./2.*mlamgl[25] + 
   mlamgl[48] + 5./16.*mlamgl[12] + 23./16. + 2*mlamgl[47];
   mlamgl[47]=mlamgl[39]*mlamgl[47];
   mlamgl[48]=51./2. - mlamgl[12];
   mlamgl[44]= - 9 + mlamgl[44];
   mlamgl[44]=mlamgl[9]*mlamgl[44];
   mlamgl[41]= - mlamgl[41] + mlamgl[44] + 1./2.*mlamgl[48] + 4*
   mlamgl[25];
   mlamgl[41]=mlamgl[41]*mlamgl[45];
   mlamgl[44]=mlamgl[21]*mlamgl[39];
   mlamgl[41]= - 2*mlamgl[44] + mlamgl[41];
   mlamgl[41]=MMt*mlamgl[41];
   mlamgl[44]=2*mlamgl[26] - mlamgl[5];
   mlamgl[44]=mlamgl[44]*mlamgl[45];
   mlamgl[41]=mlamgl[41] + 2*mlamgl[44] + mlamgl[47];
   mlamgl[41]=mlamgl[2]*mlamgl[41];
   mlamgl[44]= - 1./2.*mlamgl[29] - 6*mlamgl[20] + mlamgl[21];
   mlamgl[44]=mlamgl[39]*mlamgl[44];
   mlamgl[41]=mlamgl[44] + mlamgl[41];
   mlamgl[41]=MMt*mlamgl[41];
   mlamgl[40]=1./2.*mlamgl[40] + mlamgl[41];
   mlamgl[40]=MMt*mlamgl[40];
   mlamgl[37]=1./2.*mlamgl[37] + mlamgl[40];
   mlamgl[37]=MMt*mlamgl[37];
   mlamgl[40]=1./2.*mlamgl[8];
   mlamgl[41]=1./4.*mlamgl[17] + mlamgl[40] - 1;
   mlamgl[38]=mlamgl[41]*mlamgl[38];
   mlamgl[41]=67./8. + mlamgl[43];
   mlamgl[41]=mlamgl[41]*mlamgl[43];
   mlamgl[36]=mlamgl[36] - 1;
   mlamgl[40]=mlamgl[36]*mlamgl[40];
   mlamgl[43]= - 51./4. + 13*mlamgl[16];
   mlamgl[43]=mlamgl[43]*mlamgl[46];
   mlamgl[36]=mlamgl[8]*mlamgl[36];
   mlamgl[36]=1 + mlamgl[36];
   mlamgl[36]=mlamgl[38] + 1./2.*mlamgl[36];
   mlamgl[36]=mlamgl[35]*mlamgl[36];
   mlamgl[44]=9./2.*mlamgl[19] + mlamgl[27];
   mlamgl[36]=3*mlamgl[44] + 1./2.*mlamgl[30] + mlamgl[36];
   mlamgl[36]=MMH*mlamgl[36];
   mlamgl[44]=mlamgl[22] + mlamgl[33];
   mlamgl[36]=mlamgl[36] + mlamgl[43] + mlamgl[38] + mlamgl[40] + 11./2.
   *mlamgl[31] + 81./2.*S2 + mlamgl[42] - 131./8. + mlamgl[41] - 9./2.*
   mlamgl[44];
   mlamgl[38]=mlamgl[39]*MMH;
   mlamgl[36]=mlamgl[36]*mlamgl[38];
   mlamgl[40]=mlamgl[5] - mlamgl[6] + 1./2.*mlamgl[4];
   mlamgl[40]=mlamgl[39]*mlamgl[40];
   mlamgl[36]=mlamgl[36] + 3*mlamgl[40];
   mlamgl[36]=MMH*mlamgl[36];
   mlamgl[36]=1./16.*mlamgl[36] + mlamgl[37];
   mlamgl[37]=MMt*pow(mlamgl[2],2);
   mlamgl[37]=7*mlamgl[2] + 3*mlamgl[37];
   mlamgl[37]=MMt*mlamgl[37];
   mlamgl[37]= - 17./8. + mlamgl[37];
   mlamgl[37]=MMt*mlamgl[39]*mlamgl[37];
   mlamgl[37]=1./4.*mlamgl[38] + mlamgl[37];
   mlamgl[37]=MMt*mlamgl[37];
   mlamgl[38]=mlamgl[39]*pow(MMH,2);
   mlamgl[37]= - 23./24.*mlamgl[38] + mlamgl[37];
   mlamgl[37]=mlamgl[37]*pow(Pi,2);
   mlamgl[36]=3*mlamgl[36] + 1./8.*mlamgl[37];

      return mlamgl[36]*pow(mlamgl[3],2);
}
