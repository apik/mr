#include <tt.hpp>
std::complex<long double> tt::mgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mttgl[46];

    mttgl[1]=pow(SW,-1);
    mttgl[2]=pow(MMH,-1);
    mttgl[3]=pow(MMW,-1);
    mttgl[4]=pow(MMt,-1);
    mttgl[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    mttgl[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    mttgl[7]=Tsil::I2(0,0,MMH,mu2);
    mttgl[8]=Tsil::I2(0,0,MMt,mu2);
    mttgl[9]=Tsil::B(MMH,MMt,MMt,mu2);
    mttgl[10]=Tsil::B(MMt,MMt,MMH,mu2);
    mttgl[11]=Tsil::A(MMH,mu2);
    mttgl[12]=Tsil::A(MMt,mu2);
    mttgl[13]=Tsil::B(0,MMH,MMt,mu2);
    mttgl[14]=Tsil::B(0,0,MMH,mu2);
    mttgl[15]=Tsil::B(0,0,MMt,mu2);
    mttgl[16]=Tsil::Beps(MMH,MMt,MMt,mu2);
    mttgl[17]=Tsil::Aeps(MMH,mu2);
    mttgl[18]=Tsil::Aeps(MMt,mu2);
    mttgl[19]=prot0ttHt->Tuxv(0);
    mttgl[20]=protHHttH->M(0);
    mttgl[21]=protHttHt->M(0);
    mttgl[22]=protHHttH->Uxzuv(0);
    mttgl[23]=protHttHt->Txuv(0);
    mttgl[24]=protHHttH->Suxv(0);
    mttgl[25]=protHHttH->Uuyxv(0);
    mttgl[26]=protWt000->Tyzv(0);
    mttgl[27]=Tfin(MMH,0,0);
    mttgl[28]=Mfin(MMH,0,MMt,MMt,0);
    mttgl[29]=Mfin(MMH,0,MMt,0,0);
    mttgl[30]=Mfin(0,MMH,MMt,MMt,0);
    mttgl[31]=Mfin(0,MMH,0,MMt,0);
    mttgl[32]=Mfin(0,0,MMt,MMt,MMH);
    mttgl[33]=Ufin(MMH,MMt,0,0);
    mttgl[34]=Ufin(MMt,MMH,0,0);
    mttgl[35]=Ufin(MMt,0,0,MMH);
   mttgl[36]=mttgl[12]*pow(mttgl[2],2);
   mttgl[37]=mttgl[15]*mttgl[2];
   mttgl[38]=mttgl[2]*mttgl[10];
   mttgl[38]= - 3./2.*mttgl[36] + mttgl[38] - 5./8.*mttgl[37];
   mttgl[38]=mttgl[12]*mttgl[38];
   mttgl[39]=3*mttgl[2];
   mttgl[40]=5./2. + mttgl[10];
   mttgl[40]=mttgl[40]*mttgl[39];
   mttgl[40]= - 3./2.*mttgl[37] + mttgl[40] + mttgl[21];
   mttgl[41]=mttgl[9]*mttgl[2];
   mttgl[42]= - 3 - mttgl[10];
   mttgl[42]=mttgl[42]*mttgl[41];
   mttgl[36]=mttgl[10]*mttgl[36];
   mttgl[36]=3./2.*mttgl[42] + 1./2.*mttgl[40] - 18*mttgl[36];
   mttgl[36]=MMt*mttgl[36];
   mttgl[40]= - 11./4. - mttgl[10];
   mttgl[40]=mttgl[40]*mttgl[39];
   mttgl[37]=mttgl[40] + 1./4.*mttgl[37];
   mttgl[37]=1./2.*mttgl[37] + mttgl[41];
   mttgl[37]=mttgl[11]*mttgl[37];
   mttgl[37]=mttgl[37] + mttgl[23];
   mttgl[40]=mttgl[12]*mttgl[2];
   mttgl[41]= - 19./2. + 5*mttgl[10];
   mttgl[41]=3*mttgl[41] + 11./4.*mttgl[15];
   mttgl[41]=3./4.*mttgl[9] + 1./8.*mttgl[41] - 9*mttgl[40];
   mttgl[41]=mttgl[9]*mttgl[41];
   mttgl[42]=1 + mttgl[15];
   mttgl[42]=mttgl[15]*mttgl[42];
   mttgl[42]=mttgl[42] - mttgl[13];
   mttgl[43]=9./4.*mttgl[6] - 169./32.*mttgl[17] - 6*mttgl[18];
   mttgl[43]=mttgl[2]*mttgl[43];
   mttgl[44]=5*mttgl[33];
   mttgl[45]=13*mttgl[16] - 57./16. - mttgl[44];
   mttgl[45]=1./4.*mttgl[45] - 21*mttgl[10];
   mttgl[36]=mttgl[36] + 1./16.*mttgl[27] - 1./4.*mttgl[22] + mttgl[41]
    + 3*mttgl[38] + 75./64.*mttgl[26] + 1./8.*mttgl[45] + mttgl[43] + 3.
   /128.*mttgl[42] + 1./2.*mttgl[37];
   mttgl[36]=MMt*mttgl[36];
   mttgl[37]=1./4.*mttgl[15] + 7./8. + mttgl[10];
   mttgl[37]= - 9./2.*mttgl[9] + 3*mttgl[37] + 79./2.*mttgl[40];
   mttgl[37]=mttgl[11]*mttgl[37];
   mttgl[38]=mttgl[9]*mttgl[12];
   mttgl[40]=3*mttgl[10];
   mttgl[41]=9*mttgl[14] + 3./2.*mttgl[15] - 251./16. + mttgl[40];
   mttgl[41]=mttgl[12]*mttgl[41];
   mttgl[37]=1./8.*mttgl[37] + 3./4.*mttgl[6] + 33./4.*mttgl[38] + 1./4.
   *mttgl[41] + 129./32.*mttgl[17] + mttgl[18];
   mttgl[36]=1./2.*mttgl[37] + mttgl[36];
   mttgl[36]=MMt*mttgl[36];
   mttgl[37]= - 13./2. - mttgl[14];
   mttgl[37]=3*mttgl[37] + mttgl[9];
   mttgl[37]=mttgl[9]*mttgl[37];
   mttgl[41]= - mttgl[27] - 1 - mttgl[13];
   mttgl[41]=MMH*mttgl[4]*mttgl[41];
   mttgl[37]=mttgl[35] + mttgl[37] + mttgl[41];
   mttgl[41]=mttgl[11] - mttgl[7] + mttgl[12] + mttgl[18];
   mttgl[42]=1./2.*mttgl[4];
   mttgl[41]=mttgl[42]*mttgl[41];
   mttgl[42]=1./2.*MMt;
   mttgl[43]= - 9*mttgl[20] + mttgl[32] + mttgl[30] + mttgl[28] - 
   mttgl[21];
   mttgl[43]=mttgl[43]*mttgl[42];
   mttgl[45]= - mttgl[23] - 7./4.*mttgl[16] - 33./8. + mttgl[25];
   mttgl[37]=mttgl[43] + 3./2.*mttgl[34] + 3./4.*mttgl[22] + 5./4.*
   mttgl[13] + 9./4.*mttgl[14] + 3*mttgl[45] + mttgl[41] + 1./2.*
   mttgl[37];
   mttgl[41]=1./4.*MMH;
   mttgl[37]=mttgl[37]*mttgl[41];
   mttgl[43]=71./2. + mttgl[44];
   mttgl[43]=35./2.*mttgl[16] + 1./2.*mttgl[43] - 9*mttgl[25];
   mttgl[43]= - 3./2.*mttgl[14] - 5*mttgl[23] + 5./2.*mttgl[26] + 1./2.
   *mttgl[43] + 9*mttgl[10];
   mttgl[40]=3./4.*mttgl[14] + 55./2. - mttgl[40];
   mttgl[40]=1./2.*mttgl[40] - mttgl[9];
   mttgl[40]=mttgl[9]*mttgl[40];
   mttgl[44]=mttgl[31] + mttgl[29];
   mttgl[44]=3*mttgl[20] - mttgl[21] + 1./4.*mttgl[44];
   mttgl[44]=MMt*mttgl[44];
   mttgl[40]=mttgl[44] - 11./8.*mttgl[34] + 3./16.*mttgl[27] - 3./32.*
   mttgl[13] + 1./4.*mttgl[43] + mttgl[40];
   mttgl[40]=MMt*mttgl[40];
   mttgl[43]=mttgl[12]*mttgl[4];
   mttgl[44]=mttgl[43] - mttgl[14];
   mttgl[44]=65./2. - 3*mttgl[44];
   mttgl[44]=mttgl[12]*mttgl[44];
   mttgl[38]= - 3./2.*mttgl[5] - 5./2.*mttgl[7] - 13*mttgl[6] - 23*
   mttgl[38] + mttgl[44] + 29./2.*mttgl[17] + 37*mttgl[18];
   mttgl[44]=33./2. + mttgl[14];
   mttgl[43]=3./2.*mttgl[44] - mttgl[43];
   mttgl[43]=1./8.*mttgl[43] + mttgl[9];
   mttgl[43]=mttgl[11]*mttgl[43];
   mttgl[37]=mttgl[37] + mttgl[40] + 1./8.*mttgl[38] + mttgl[43];
   mttgl[37]=mttgl[37]*mttgl[41];
   mttgl[38]=pow(MMt,2);
   mttgl[40]=35*MMt - 23./4.*MMH;
   mttgl[40]=MMH*mttgl[40];
   mttgl[38]= - 115./2.*mttgl[38] + mttgl[40];
   mttgl[38]=mttgl[19]*mttgl[38];
   mttgl[40]= - 33*mttgl[12] - 29./8.*mttgl[11];
   mttgl[40]=mttgl[11]*mttgl[40];
   mttgl[38]=mttgl[40] + mttgl[38];
   mttgl[39]=MMt*mttgl[39];
   mttgl[39]=5./16. + mttgl[39];
   mttgl[39]=mttgl[8]*mttgl[39]*mttgl[42];
   mttgl[40]=pow(mttgl[12],2);
   mttgl[41]= - 5./4.*MMt - MMH;
   mttgl[41]=mttgl[24]*mttgl[41];
   mttgl[36]=mttgl[39] + 3./16.*mttgl[41] + mttgl[37] + mttgl[36] + 
   mttgl[40] + 1./16.*mttgl[38];

      return mttgl[36]*pow(mttgl[3],2)*pow(mttgl[1],4);
}
