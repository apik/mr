#include <HH.hpp>
std::complex<long double> HH::mgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHHgl[47];

    mHHgl[1]=pow(SW,-1);
    mHHgl[2]=pow(MMH,-1);
    mHHgl[3]=pow(MMW,-1);
    mHHgl[4]=Tsil::I2(MMH,MMH,MMH,mu2);
    mHHgl[5]=Tsil::I2(MMH,MMt,MMt,mu2);
    mHHgl[6]=Tsil::B(MMH,MMH,MMH,mu2);
    mHHgl[7]=Tsil::B(MMt,MMt,MMH,mu2);
    mHHgl[8]=log(pow(mu2,-1)*MMt);
    mHHgl[9]=log(pow(mu2,-1)*MMH);
    mHHgl[10]=Tsil::B(MMH,MMt,MMt,mu2);
    mHHgl[11]=Tsil::B(0,MMt,MMH,mu2);
    mHHgl[12]=Tsil::B(0,0,MMt,mu2);
    mHHgl[13]=Tsil::Beps(MMH,MMH,MMH,mu2);
    mHHgl[14]=Tsil::Beps(MMt,MMt,MMH,mu2);
    mHHgl[15]=log(MMH);
    mHHgl[16]=log(MMt);
    mHHgl[17]=prottttt0->Suxv(0);
    mHHgl[18]=protHHHHH->M(0);
    mHHgl[19]=protHtHtt->M(0);
    mHHgl[20]=protttttH->M(0);
    mHHgl[21]=protHHHHH->Uzxyv(0);
    mHHgl[22]=protHtHtt->Uzxyv(0);
    mHHgl[23]=protHtHtt->Uuyxv(0);
    mHHgl[24]=protHtHtt->Tyzv(0);
    mHHgl[25]=protHtHtt->Svyz(0);
    mHHgl[26]=Mfin(0,MMH,0,MMH,0);
    mHHgl[27]=Mfin(0,MMt,0,MMt,MMt);
    mHHgl[28]=Mfin(0,MMt,0,MMt,0);
    mHHgl[29]=Mfin(0,0,0,0,MMH);
    mHHgl[30]=Tfin(MMH,0,0);
    mHHgl[31]=Tfin(MMt,0,0);
    mHHgl[32]=Ufin(MMH,MMH,0,0);
    mHHgl[33]=Ufin(MMt,MMt,0,0);
    mHHgl[34]=1/(4*MMt - MMH);
   mHHgl[35]=5*mHHgl[8];
   mHHgl[36]= - 25./16. - mHHgl[10];
   mHHgl[36]=mHHgl[36]*mHHgl[35];
   mHHgl[37]=5*mHHgl[10];
   mHHgl[38]=1./4.*mHHgl[8];
   mHHgl[39]= - 7./8.*mHHgl[7] + mHHgl[38] + 3./2. - mHHgl[37];
   mHHgl[39]=mHHgl[7]*mHHgl[39];
   mHHgl[40]= - mHHgl[5] + 2*mHHgl[25];
   mHHgl[40]=mHHgl[2]*mHHgl[40];
   mHHgl[40]=mHHgl[40] + mHHgl[23];
   mHHgl[41]=pow(Pi,2);
   mHHgl[42]=mHHgl[16]*mHHgl[8];
   mHHgl[43]= - 2 - mHHgl[7];
   mHHgl[43]=2*mHHgl[43] + mHHgl[15];
   mHHgl[43]=mHHgl[9]*mHHgl[43];
   mHHgl[44]=mHHgl[7] + mHHgl[8];
   mHHgl[44]=1 - 5./4.*mHHgl[44];
   mHHgl[44]=mHHgl[12]*mHHgl[44];
   mHHgl[36]=mHHgl[43] + 19./32.*mHHgl[41] - 3./8.*mHHgl[31] + 5./16.*
   mHHgl[11] + 5./4.*mHHgl[33] + 43./8.*mHHgl[42] + mHHgl[39] + 
   mHHgl[36] - 5./2.*mHHgl[24] + 23./2. + 4*mHHgl[10] + 2*mHHgl[40] - 
   13./4.*mHHgl[14] + mHHgl[44];
   mHHgl[36]=mHHgl[2]*mHHgl[36];
   mHHgl[39]=1./2.*mHHgl[31];
   mHHgl[40]=1./8.*mHHgl[41] - mHHgl[39] - 1./2.*mHHgl[11];
   mHHgl[40]=mHHgl[40]*pow(mHHgl[2],2);
   mHHgl[43]=1./2.*mHHgl[42] - 9*mHHgl[8] + 51./4. + 4*mHHgl[24];
   mHHgl[43]=mHHgl[2]*mHHgl[43];
   mHHgl[43]= - 2*mHHgl[20] + mHHgl[43];
   mHHgl[43]=mHHgl[2]*mHHgl[43];
   mHHgl[40]=mHHgl[43] + mHHgl[40];
   mHHgl[40]=MMt*mHHgl[40];
   mHHgl[36]=mHHgl[40] - 1./2.*mHHgl[28] + mHHgl[20] - 6*mHHgl[19] + 
   mHHgl[36];
   mHHgl[36]=MMt*mHHgl[36];
   mHHgl[40]= - 5./2.*mHHgl[25] + 11*mHHgl[5] - 1./2.*mHHgl[17];
   mHHgl[40]=mHHgl[2]*mHHgl[40];
   mHHgl[43]= - 1 - 3*mHHgl[7];
   mHHgl[43]=17./4.*mHHgl[15] + 13./2.*mHHgl[43] - 11*mHHgl[16];
   mHHgl[43]=mHHgl[9]*mHHgl[43];
   mHHgl[44]= - mHHgl[27] + 1./2.*mHHgl[20] + 9*mHHgl[19];
   mHHgl[44]=MMH*mHHgl[44];
   mHHgl[40]=mHHgl[44] + mHHgl[40] + mHHgl[43];
   mHHgl[37]=133./4. + mHHgl[37];
   mHHgl[37]=mHHgl[37]*mHHgl[38];
   mHHgl[38]=1./4.*mHHgl[7];
   mHHgl[43]= - mHHgl[38] + 3./2.*mHHgl[8] + 81./4. + 7*mHHgl[10];
   mHHgl[43]=mHHgl[43]*mHHgl[38];
   mHHgl[44]= - 7./2.*mHHgl[7] + 1 - 5./2.*mHHgl[8];
   mHHgl[44]=mHHgl[6]*mHHgl[44];
   mHHgl[44]= - mHHgl[13] + mHHgl[44] + mHHgl[22];
   mHHgl[45]=mHHgl[33] + mHHgl[42] - mHHgl[14];
   mHHgl[46]=mHHgl[12]*mHHgl[7];
   mHHgl[36]=mHHgl[36] + 1./8.*mHHgl[46] + 3./32.*mHHgl[41] + 9./16.*
   mHHgl[11] + mHHgl[43] + mHHgl[37] - 3./2.*mHHgl[24] - 349./64. - 
   mHHgl[10] + 1./4.*mHHgl[40] - 5./8.*mHHgl[45] + 9./4.*mHHgl[44];
   mHHgl[37]=pow(mHHgl[1],4);
   mHHgl[40]=pow(mHHgl[3],2);
   mHHgl[43]=mHHgl[37]*mHHgl[40];
   mHHgl[36]=MMt*mHHgl[43]*mHHgl[36];
   mHHgl[38]=mHHgl[38] + 1./2.*mHHgl[8];
   mHHgl[44]=3./2. - mHHgl[10] - mHHgl[38];
   mHHgl[44]=mHHgl[7]*mHHgl[44];
   mHHgl[35]=mHHgl[44] - 47./4. + mHHgl[35];
   mHHgl[44]=mHHgl[11] + mHHgl[23];
   mHHgl[35]=mHHgl[39] - 3./16.*mHHgl[42] + 1./4.*mHHgl[35] + 3*
   mHHgl[13] - 3./4.*mHHgl[44];
   mHHgl[39]=3./2.*mHHgl[9];
   mHHgl[44]= - 3./2.*mHHgl[15] - 9./4.*mHHgl[16] + 61./8. + mHHgl[7];
   mHHgl[44]=mHHgl[44]*mHHgl[39];
   mHHgl[45]=9*mHHgl[6];
   mHHgl[46]=7./2. - mHHgl[8];
   mHHgl[46]=3./8.*mHHgl[46] + mHHgl[7];
   mHHgl[46]=mHHgl[46]*mHHgl[45];
   mHHgl[35]=mHHgl[46] + 9./4.*mHHgl[14] - 9*mHHgl[22] + mHHgl[44] + 3*
   mHHgl[35] - 1./2.*mHHgl[41];
   mHHgl[35]=MMH*mHHgl[40]*mHHgl[35];
   mHHgl[44]=3*mHHgl[40];
   mHHgl[46]= - 3./2.*mHHgl[5] - mHHgl[17];
   mHHgl[46]=mHHgl[46]*mHHgl[44];
   mHHgl[35]=mHHgl[46] + mHHgl[35];
   mHHgl[35]=mHHgl[37]*mHHgl[35];
   mHHgl[35]=1./2.*mHHgl[35] + 3*mHHgl[36];
   mHHgl[35]=MMt*mHHgl[35];
   mHHgl[36]=1./8.*mHHgl[42];
   mHHgl[38]=mHHgl[38] - 1;
   mHHgl[38]=mHHgl[38]*mHHgl[7];
   mHHgl[42]=mHHgl[38] + 11*mHHgl[30] - 53./4. - mHHgl[8];
   mHHgl[42]= - 9./2.*mHHgl[21] + mHHgl[36] + 1./2.*mHHgl[42] + 9*
   mHHgl[13];
   mHHgl[46]= - 39./4. + 7*mHHgl[15];
   mHHgl[39]=mHHgl[46]*mHHgl[39];
   mHHgl[46]=53./4. + 3*mHHgl[9];
   mHHgl[46]=1./2.*mHHgl[46] + 3*mHHgl[6];
   mHHgl[45]=mHHgl[45]*mHHgl[46];
   mHHgl[39]=mHHgl[45] - 27./2.*mHHgl[32] + mHHgl[39] + 3*mHHgl[42] - 5.
   /4.*mHHgl[41];
   mHHgl[39]=mHHgl[40]*mHHgl[39];
   mHHgl[38]=mHHgl[38] + 1 - mHHgl[8];
   mHHgl[38]=mHHgl[34]*mHHgl[38];
   mHHgl[38]=mHHgl[38] + mHHgl[29];
   mHHgl[36]=mHHgl[34]*mHHgl[36];
   mHHgl[36]=mHHgl[36] + 27./2.*mHHgl[18] + 3*mHHgl[26] + 1./2.*
   mHHgl[38];
   mHHgl[36]=MMH*mHHgl[36]*mHHgl[44];
   mHHgl[36]=mHHgl[36] + mHHgl[39];
   mHHgl[36]=MMH*mHHgl[37]*mHHgl[36];
   mHHgl[37]=mHHgl[4]*mHHgl[43];
   mHHgl[36]=9./2.*mHHgl[37] + mHHgl[36];
   mHHgl[36]=MMH*mHHgl[36];

      return mHHgl[35] + 1./16.*mHHgl[36];
}
