#include <HH.hpp>
std::complex<long double> HH::mgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mHHgl[49];

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
    mHHgl[26]=prot0H0H0->M(0);
    mHHgl[27]=prot0t0tt->M(0);
    mHHgl[28]=prot0t0t0->M(0);
    mHHgl[29]=prot0000H->M(0);
    mHHgl[30]=prot0H0H0->Uuyxv(0);
    mHHgl[31]=prot0t0t0->Uuyxv(0);
    mHHgl[32]=prot0H0H0->Tyzv(0);
    mHHgl[33]=prot0t0t0->Tyzv(0);
    mHHgl[34]=1/(4*MMt - MMH);
   mHHgl[35]=7./4.*mHHgl[7];
   mHHgl[36]=1./2.*mHHgl[8];
   mHHgl[37]= - mHHgl[35] + 3 + mHHgl[36];
   mHHgl[37]=mHHgl[7]*mHHgl[37];
   mHHgl[38]=pow(Pi,2);
   mHHgl[39]=mHHgl[16]*mHHgl[8];
   mHHgl[37]= - 13./2.*mHHgl[14] + 19./16.*mHHgl[38] + mHHgl[37] + 43./
   4.*mHHgl[39] - 125./8.*mHHgl[8] + 23 - 5*mHHgl[24];
   mHHgl[40]=2*mHHgl[2];
   mHHgl[41]= - mHHgl[5] + 2*mHHgl[25];
   mHHgl[41]=mHHgl[41]*mHHgl[40];
   mHHgl[42]= - 2 - mHHgl[7];
   mHHgl[42]=mHHgl[9]*mHHgl[42];
   mHHgl[42]=mHHgl[42] + mHHgl[23];
   mHHgl[43]=mHHgl[15]*mHHgl[9];
   mHHgl[44]=mHHgl[7] + mHHgl[8];
   mHHgl[44]=4 - 5*mHHgl[44];
   mHHgl[44]=mHHgl[10]*mHHgl[44];
   mHHgl[45]=5./4.*mHHgl[8];
   mHHgl[46]=mHHgl[45] - 1;
   mHHgl[47]= - 5./4.*mHHgl[7] - mHHgl[46];
   mHHgl[47]=mHHgl[12]*mHHgl[47];
   mHHgl[37]=5./4.*mHHgl[31] + mHHgl[43] + mHHgl[41] + 5./16.*mHHgl[11]
    + mHHgl[47] - 3./8.*mHHgl[33] + 1./2.*mHHgl[37] + mHHgl[44] + 2*
   mHHgl[42];
   mHHgl[37]=mHHgl[2]*mHHgl[37];
   mHHgl[41]= - mHHgl[33] + mHHgl[39] - mHHgl[11];
   mHHgl[41]=1./8.*mHHgl[38] - 9*mHHgl[8] + 51./4. + 4*mHHgl[24] + 1./2.
   *mHHgl[41];
   mHHgl[41]=mHHgl[41]*pow(mHHgl[2],2);
   mHHgl[40]= - mHHgl[20]*mHHgl[40];
   mHHgl[40]=mHHgl[41] + mHHgl[40];
   mHHgl[40]=MMt*mHHgl[40];
   mHHgl[37]=mHHgl[40] + mHHgl[20] - 6*mHHgl[19] - 1./2.*mHHgl[28] + 
   mHHgl[37];
   mHHgl[37]=MMt*mHHgl[37];
   mHHgl[40]=1./2.*mHHgl[7];
   mHHgl[41]= - 21*mHHgl[6] + 27./2. + mHHgl[8];
   mHHgl[41]=3*mHHgl[41] - mHHgl[40];
   mHHgl[42]=1./4.*mHHgl[7];
   mHHgl[41]=mHHgl[41]*mHHgl[42];
   mHHgl[44]=1 - 5./2.*mHHgl[8];
   mHHgl[44]=mHHgl[6]*mHHgl[44];
   mHHgl[44]=mHHgl[44] + mHHgl[22];
   mHHgl[47]=mHHgl[39] - mHHgl[14];
   mHHgl[41]=3./16.*mHHgl[38] + mHHgl[41] + 133./8.*mHHgl[8] - 349./32.
    - 3*mHHgl[24] - 5./4.*mHHgl[47] + 9./2.*mHHgl[44];
   mHHgl[35]=mHHgl[35] + mHHgl[46];
   mHHgl[35]=mHHgl[10]*mHHgl[35];
   mHHgl[44]=9*mHHgl[19] - mHHgl[27];
   mHHgl[44]=1./8.*mHHgl[20] + 1./4.*mHHgl[44];
   mHHgl[44]=MMH*mHHgl[44];
   mHHgl[46]=11*mHHgl[5] - 5./2.*mHHgl[25];
   mHHgl[46]= - 1./8.*mHHgl[17] + 1./4.*mHHgl[46];
   mHHgl[46]=mHHgl[2]*mHHgl[46];
   mHHgl[47]=mHHgl[12]*mHHgl[7];
   mHHgl[48]= - 39./2.*mHHgl[7] - 13./2. - 11*mHHgl[16];
   mHHgl[48]=mHHgl[9]*mHHgl[48];
   mHHgl[35]=mHHgl[37] - 5./8.*mHHgl[31] + 17./16.*mHHgl[43] + 1./4.*
   mHHgl[48] + 9./16.*mHHgl[11] + 1./8.*mHHgl[47] - 9./4.*mHHgl[13] + 1.
   /2.*mHHgl[41] + mHHgl[35] + mHHgl[46] + mHHgl[44];
   mHHgl[37]=pow(mHHgl[3],2)*pow(mHHgl[1],4);
   mHHgl[35]=MMt*mHHgl[37]*mHHgl[35];
   mHHgl[41]=3*mHHgl[6];
   mHHgl[44]=3 - mHHgl[8];
   mHHgl[44]= - 1./16.*mHHgl[7] + 1./8.*mHHgl[44] + mHHgl[41];
   mHHgl[44]=mHHgl[7]*mHHgl[44];
   mHHgl[46]=7./2. - mHHgl[8];
   mHHgl[46]=mHHgl[6]*mHHgl[46];
   mHHgl[44]=mHHgl[44] - 3./16.*mHHgl[39] + 9./8.*mHHgl[46] + mHHgl[45]
    - 47./16. - 3*mHHgl[22];
   mHHgl[45]=mHHgl[43] + mHHgl[11] + mHHgl[23] - mHHgl[14];
   mHHgl[46]=61./2. - 9*mHHgl[16];
   mHHgl[46]=1./4.*mHHgl[46] + mHHgl[7];
   mHHgl[46]=mHHgl[9]*mHHgl[46];
   mHHgl[46]=mHHgl[46] + mHHgl[33];
   mHHgl[47]=mHHgl[10]*mHHgl[7];
   mHHgl[44]=9*mHHgl[13] - 3./4.*mHHgl[47] + 3*mHHgl[44] - 1./2.*
   mHHgl[38] + 3./2.*mHHgl[46] - 9./4.*mHHgl[45];
   mHHgl[44]=MMH*mHHgl[44];
   mHHgl[44]= - 3*mHHgl[17] - 9./2.*mHHgl[5] + mHHgl[44];
   mHHgl[44]=mHHgl[44]*mHHgl[37];
   mHHgl[35]=1./2.*mHHgl[44] + 3*mHHgl[35];
   mHHgl[35]=MMt*mHHgl[35];
   mHHgl[44]=53./8. + mHHgl[41];
   mHHgl[44]=mHHgl[44]*mHHgl[41];
   mHHgl[36]=mHHgl[42] + mHHgl[36] - 1;
   mHHgl[40]=mHHgl[36]*mHHgl[40];
   mHHgl[42]= - 9*mHHgl[21] - 53./4. - mHHgl[8];
   mHHgl[40]=mHHgl[40] + 1./8.*mHHgl[39] + 1./2.*mHHgl[42] + mHHgl[44];
   mHHgl[41]= - 13./4. + mHHgl[41];
   mHHgl[41]=mHHgl[9]*mHHgl[41];
   mHHgl[38]=21./2.*mHHgl[43] + 9./2.*mHHgl[41] - 27./2.*mHHgl[30] + 27
   *mHHgl[13] + 33./2.*mHHgl[32] + 3*mHHgl[40] - 5./4.*mHHgl[38];
   mHHgl[38]=mHHgl[38]*pow(MMH,2);
   mHHgl[36]=mHHgl[7]*mHHgl[36];
   mHHgl[36]=mHHgl[36] + 1./4.*mHHgl[39] + 1 - mHHgl[8];
   mHHgl[36]=mHHgl[34]*mHHgl[36];
   mHHgl[36]=mHHgl[36] + mHHgl[29];
   mHHgl[36]=9*mHHgl[26] + 81./2.*mHHgl[18] + 3./2.*mHHgl[36];
   mHHgl[36]=mHHgl[36]*pow(MMH,3);
   mHHgl[39]=mHHgl[4]*MMH;
   mHHgl[36]=9./2.*mHHgl[39] + mHHgl[38] + mHHgl[36];
   mHHgl[36]=mHHgl[36]*mHHgl[37];

      return mHHgl[35] + 1./16.*mHHgl[36];
}
