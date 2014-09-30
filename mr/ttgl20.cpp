#include <tt.hpp>
std::complex<long double> tt::mgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mttgl[51];

    mttgl[1]=pow(SW,-1);
    mttgl[2]=pow(MMH,-1);
    mttgl[3]=pow(MMW,-1);
    mttgl[4]=pow(MMt,-1);
    mttgl[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    mttgl[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    mttgl[7]=Tsil::I2(0,0,MMH,mu2);
    mttgl[8]=Tsil::I2(0,0,MMt,mu2);
    mttgl[9]=Tsil::B(MMH,MMH,MMH,mu2);
    mttgl[10]=Tsil::B(MMH,MMt,MMt,mu2);
    mttgl[11]=Tsil::A(MMH,mu2);
    mttgl[12]=Tsil::A(MMt,mu2);
    mttgl[13]=Tsil::B(MMt,MMt,MMH,mu2);
    mttgl[14]=std::real(Tsil::B(0,0,MMH,mu2));
    mttgl[15]=std::real(Tsil::B(0,0,MMt,mu2));
    mttgl[16]=Tsil::B(0,MMH,MMt,mu2);
    mttgl[17]=Tsil::Beps(MMH,MMt,MMt,mu2);
    mttgl[18]=Tsil::Aeps(MMH,mu2);
    mttgl[19]=Tsil::Aeps(MMt,mu2);
    mttgl[20]=prot0ttHt->Tuxv(0);
    mttgl[21]=protHHttH->M(0);
    mttgl[22]=protHttHt->M(0);
    mttgl[23]=protHHttH->Uxzuv(0);
    mttgl[24]=protHttHt->Txuv(0);
    mttgl[25]=protHHttH->Suxv(0);
    mttgl[26]=protHHttH->Uuyxv(0);
    mttgl[27]=protWt000->Tyzv(0);
    mttgl[28]=protH0tt0->M(0);
    mttgl[29]=protH0t00->M(0);
    mttgl[30]=prot0Htt0->M(0);
    mttgl[31]=prot0H0t0->M(0);
    mttgl[32]=prot00ttH->M(0);
    mttgl[33]=protH0t00->Uxzuv(0);
    mttgl[34]=protH0t00->Uzxyv(0);
    mttgl[35]=protH0tt0->Uuyxv(0);
    mttgl[36]=protH0t00->Txuv(0);
   mttgl[37]=3*mttgl[10];
   mttgl[38]= - 19./4. + mttgl[10];
   mttgl[38]=mttgl[38]*mttgl[37];
   mttgl[39]=3./32.*mttgl[16];
   mttgl[40]=5*mttgl[33];
   mttgl[41]= - 57./16. - mttgl[40];
   mttgl[38]= - mttgl[23] + 13./8.*mttgl[17] - mttgl[39] + 1./8.*
   mttgl[41] + mttgl[38];
   mttgl[41]=1./2.*mttgl[11];
   mttgl[42]= - 33./8. + mttgl[10];
   mttgl[42]=mttgl[42]*mttgl[41];
   mttgl[43]=mttgl[10]*mttgl[12];
   mttgl[44]=pow(mttgl[12],2);
   mttgl[45]=mttgl[2]*mttgl[44];
   mttgl[42]= - 9./2.*mttgl[45] + 9./4.*mttgl[6] - 6*mttgl[19] - 169./
   32.*mttgl[18] - 9*mttgl[43] + mttgl[42];
   mttgl[42]=mttgl[2]*mttgl[42];
   mttgl[45]=3*mttgl[2];
   mttgl[37]=5./2. - mttgl[37];
   mttgl[37]=mttgl[37]*mttgl[45];
   mttgl[37]=mttgl[22] + mttgl[37];
   mttgl[46]=mttgl[10] - 1;
   mttgl[47]=mttgl[2]*mttgl[12];
   mttgl[47]= - 1./2.*mttgl[46] - 6*mttgl[47];
   mttgl[47]=mttgl[13]*mttgl[47]*mttgl[45];
   mttgl[48]=mttgl[15]*mttgl[2];
   mttgl[37]= - 3./4.*mttgl[48] + 1./2.*mttgl[37] + mttgl[47];
   mttgl[37]=MMt*mttgl[37];
   mttgl[47]= - 7 + 5*mttgl[10];
   mttgl[48]=1./4.*mttgl[11];
   mttgl[49]=mttgl[12] - mttgl[48];
   mttgl[49]=mttgl[2]*mttgl[49];
   mttgl[47]=1./8.*mttgl[47] + mttgl[49];
   mttgl[47]=mttgl[13]*mttgl[47];
   mttgl[49]=3./4. + 11*mttgl[10];
   mttgl[50]= - 15*mttgl[12] + mttgl[41];
   mttgl[50]=mttgl[2]*mttgl[50];
   mttgl[49]=3./16.*mttgl[15] + 1./4.*mttgl[49] + mttgl[50];
   mttgl[49]=mttgl[15]*mttgl[49];
   mttgl[37]=mttgl[37] + 1./8.*mttgl[49] + 3*mttgl[47] + 1./2.*
   mttgl[24] + 75./64.*mttgl[27] + mttgl[42] + 1./4.*mttgl[38];
   mttgl[37]=MMt*mttgl[37];
   mttgl[38]=mttgl[41] + mttgl[12];
   mttgl[41]=mttgl[13]*mttgl[38];
   mttgl[41]=mttgl[41] + mttgl[6];
   mttgl[42]=mttgl[2]*mttgl[11];
   mttgl[42]=187./16.*mttgl[42] + 9./4.*mttgl[14] + 27./4.*mttgl[9];
   mttgl[42]=mttgl[12]*mttgl[42];
   mttgl[47]=129./8.*mttgl[18] - 15./16.*mttgl[11] - 467./16.*mttgl[12]
    + 33*mttgl[43];
   mttgl[48]=mttgl[12] + mttgl[48];
   mttgl[48]=mttgl[15]*mttgl[48];
   mttgl[41]=3./8.*mttgl[48] - 15./32.*mttgl[25] + 1./4.*mttgl[47] + 
   mttgl[19] + mttgl[42] + 3./4.*mttgl[41];
   mttgl[37]=1./2.*mttgl[41] + mttgl[37];
   mttgl[37]=MMt*mttgl[37];
   mttgl[41]=1./2.*MMt;
   mttgl[42]=MMt*mttgl[45];
   mttgl[42]=5./16. + mttgl[42];
   mttgl[42]=mttgl[8]*mttgl[42]*mttgl[41];
   mttgl[45]=-pow(Pi,2);
   mttgl[45]=1./16.*mttgl[36] - 3./128.*mttgl[45] - 115./32.*mttgl[20];
   mttgl[45]=mttgl[45]*pow(MMt,2);
   mttgl[47]= - 57*mttgl[12] - 11./4.*mttgl[11];
   mttgl[47]=mttgl[11]*mttgl[47];
   mttgl[37]=mttgl[42] + mttgl[37] + mttgl[44] + 1./32.*mttgl[47] + 
   mttgl[45];
   mttgl[42]=pow(mttgl[3],2)*pow(mttgl[1],4);
   mttgl[37]=mttgl[37]*mttgl[42];
   mttgl[45]=3./8.*mttgl[14] + 9./8.*mttgl[9];
   mttgl[45]=mttgl[46]*mttgl[45];
   mttgl[40]=89./2. + mttgl[40];
   mttgl[46]=211./16. - mttgl[10];
   mttgl[46]=mttgl[10]*mttgl[46];
   mttgl[47]=mttgl[10] - 3./2.;
   mttgl[48]=mttgl[13]*mttgl[47];
   mttgl[49]=mttgl[29] + mttgl[31];
   mttgl[49]= - mttgl[22] + 3*mttgl[21] + 1./4.*mttgl[49];
   mttgl[49]=MMt*mttgl[49];
   mttgl[39]=35./4.*mttgl[20] + 3./16.*mttgl[36] - 9./8.*mttgl[26] + 
   mttgl[49] - 3./2.*mttgl[48] - 5./4.*mttgl[24] + 5./8.*mttgl[27] - 11.
   /8.*mttgl[34] + 35./16.*mttgl[17] - mttgl[39] + 1./16.*mttgl[40] + 
   mttgl[46] + mttgl[45];
   mttgl[39]=MMt*mttgl[39];
   mttgl[40]=3./4.*mttgl[14] + 9./4.*mttgl[9];
   mttgl[38]=mttgl[38]*mttgl[40];
   mttgl[40]= - mttgl[11]*mttgl[12];
   mttgl[40]= - 3*mttgl[44] + mttgl[40];
   mttgl[40]=mttgl[4]*mttgl[40];
   mttgl[44]=9 - 1./4.*mttgl[10];
   mttgl[44]=mttgl[11]*mttgl[44];
   mttgl[38]= - 3./8.*mttgl[5] + 1./4.*mttgl[40] - 3./2.*mttgl[25] - 5./
   8.*mttgl[7] - 13./4.*mttgl[6] + 37./4.*mttgl[19] + 29./8.*mttgl[18]
    + mttgl[44] + 7*mttgl[12] - 23./4.*mttgl[43] + mttgl[38];
   mttgl[40]= - 3*mttgl[14] - 9*mttgl[9];
   mttgl[40]=mttgl[47]*mttgl[40];
   mttgl[43]= - 15 + mttgl[10];
   mttgl[43]=mttgl[10]*mttgl[43];
   mttgl[44]= - mttgl[7] + mttgl[19] + mttgl[12] + mttgl[11];
   mttgl[44]=mttgl[4]*mttgl[44];
   mttgl[45]= - mttgl[36] - 1 - mttgl[16];
   mttgl[45]=MMH*mttgl[4]*mttgl[45];
   mttgl[40]=mttgl[45] + mttgl[44] + 3*mttgl[34] + mttgl[35] - 21./2.*
   mttgl[17] + 5./2.*mttgl[16] - 63./2. + mttgl[43] + mttgl[40];
   mttgl[43]=mttgl[30] + mttgl[32] - mttgl[22] + mttgl[28] - 9*
   mttgl[21];
   mttgl[41]=mttgl[41]*mttgl[43];
   mttgl[43]=mttgl[26] - mttgl[24];
   mttgl[40]= - 23./4.*mttgl[20] + 3./4.*mttgl[23] + 1./2.*mttgl[40] + 
   3*mttgl[43] + mttgl[41];
   mttgl[41]=1./4.*MMH;
   mttgl[40]=mttgl[41]*mttgl[40];
   mttgl[38]=mttgl[40] + 1./2.*mttgl[38] + mttgl[39];
   mttgl[38]=mttgl[41]*mttgl[42]*mttgl[38];

      return mttgl[37] + mttgl[38];
}
