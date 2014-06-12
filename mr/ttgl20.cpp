#include <tt.hpp>
std::complex<long double> tt::mgl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mttgl[47];

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
    mttgl[27]=protH0tt0->M(0);
    mttgl[28]=protH0t00->M(0);
    mttgl[29]=prot0Htt0->M(0);
    mttgl[30]=prot0H0t0->M(0);
    mttgl[31]=prot00ttH->M(0);
    mttgl[32]=protH0t00->Uxzuv(0);
    mttgl[33]=protH0t00->Uzxyv(0);
    mttgl[34]=protH0tt0->Uuyxv(0);
    mttgl[35]=protH0t00->Txuv(0);
   mttgl[36]=1./2.*mttgl[9];
   mttgl[37]= - 19./2. + 5*mttgl[10];
   mttgl[37]=3*mttgl[37] + 11./4.*mttgl[15];
   mttgl[37]=1./2.*mttgl[37] + 3*mttgl[9];
   mttgl[37]=mttgl[37]*mttgl[36];
   mttgl[38]=1 + mttgl[15];
   mttgl[38]=mttgl[15]*mttgl[38];
   mttgl[38]=mttgl[38] - mttgl[13];
   mttgl[39]=5./16.*mttgl[32];
   mttgl[40]=MMt*mttgl[21];
   mttgl[37]=mttgl[40] - mttgl[39] + mttgl[37] + 1./8.*mttgl[35] - 21./
   4.*mttgl[10] - 1./2.*mttgl[22] - 115./16.*mttgl[19] - 57./256. + 
   mttgl[23] + 3./64.*mttgl[38];
   mttgl[37]=MMt*mttgl[37];
   mttgl[38]=mttgl[14]*mttgl[12];
   mttgl[40]=mttgl[15]*mttgl[12];
   mttgl[41]=mttgl[10]*mttgl[12];
   mttgl[42]=mttgl[9]*mttgl[12];
   mttgl[43]= - 251./16.*mttgl[12] + 3*mttgl[6] + 5./4.*mttgl[8];
   mttgl[37]=mttgl[37] + 33./4.*mttgl[42] + 3./8.*mttgl[40] + 3./4.*
   mttgl[41] + 9./4.*mttgl[38] + 129./32.*mttgl[17] + 1./4.*mttgl[43]
    + mttgl[18];
   mttgl[43]=1./2.*MMt;
   mttgl[37]=mttgl[37]*mttgl[43];
   mttgl[44]=pow(mttgl[12],2);
   mttgl[45]=MMt*mttgl[41];
   mttgl[45]= - 1./2.*mttgl[44] - 2*mttgl[45];
   mttgl[45]=mttgl[2]*mttgl[45];
   mttgl[45]=mttgl[42] - mttgl[45];
   mttgl[46]=3./2.*mttgl[6] + mttgl[8];
   mttgl[41]=mttgl[41] + 1./2.*mttgl[46] - 2*mttgl[18];
   mttgl[46]= - 3 - mttgl[10];
   mttgl[46]=mttgl[9]*mttgl[46];
   mttgl[46]=mttgl[46] - 1./2.*mttgl[15] + 5./2. + mttgl[10];
   mttgl[46]=MMt*mttgl[46];
   mttgl[40]=3./2.*mttgl[46] - 15./8.*mttgl[40] - 169./32.*mttgl[17] + 
   3*mttgl[41] - 9*mttgl[45];
   mttgl[41]=pow(MMt,2);
   mttgl[40]=mttgl[41]*mttgl[40];
   mttgl[45]=1./4.*mttgl[15];
   mttgl[46]= - 11./4. - mttgl[10];
   mttgl[46]=3*mttgl[46] + mttgl[45];
   mttgl[46]=1./2.*mttgl[46] + mttgl[9];
   mttgl[46]=MMt*mttgl[46];
   mttgl[46]=79./16.*mttgl[12] + mttgl[46];
   mttgl[46]=mttgl[11]*mttgl[46]*mttgl[43];
   mttgl[40]=mttgl[46] + mttgl[40];
   mttgl[40]=mttgl[2]*mttgl[40];
   mttgl[46]=75./64.*mttgl[26] + 13./32.*mttgl[16];
   mttgl[41]=mttgl[41]*mttgl[46];
   mttgl[46]=mttgl[24]*MMt;
   mttgl[45]= - 3./2.*mttgl[9] + mttgl[45] + 7./8. + mttgl[10];
   mttgl[45]=MMt*mttgl[45];
   mttgl[45]= - 11*mttgl[12] + mttgl[45];
   mttgl[45]=3*mttgl[45] - 29./8.*mttgl[11];
   mttgl[45]=mttgl[11]*mttgl[45];
   mttgl[37]=mttgl[40] + 1./16.*mttgl[45] - 15./64.*mttgl[46] + 
   mttgl[44] + mttgl[37] + mttgl[41];
   mttgl[40]=pow(mttgl[3],2)*pow(mttgl[1],4);
   mttgl[37]=mttgl[40]*mttgl[37];
   mttgl[41]=3./2.*mttgl[14];
   mttgl[45]= - 11./2.*mttgl[33] + 3./4.*mttgl[35] + 9*mttgl[10] - 
   mttgl[41] + 35*mttgl[19] + 71./8. - 5*mttgl[23];
   mttgl[41]=55 + mttgl[41];
   mttgl[41]=1./2.*mttgl[41] - 3*mttgl[10];
   mttgl[41]=1./2.*mttgl[41] - mttgl[9];
   mttgl[41]=mttgl[9]*mttgl[41];
   mttgl[46]=mttgl[28] + mttgl[30];
   mttgl[46]= - mttgl[21] + 3*mttgl[20] + 1./4.*mttgl[46];
   mttgl[46]=MMt*mttgl[46];
   mttgl[39]=5./8.*mttgl[26] + 35./16.*mttgl[16] + mttgl[46] - 3./32.*
   mttgl[13] - 9./8.*mttgl[25] + mttgl[39] + 1./4.*mttgl[45] + 
   mttgl[41];
   mttgl[39]=MMt*mttgl[39];
   mttgl[41]= - mttgl[13] - 1 - mttgl[35];
   mttgl[41]=MMH*mttgl[41];
   mttgl[41]= - mttgl[7] + mttgl[12] + mttgl[18] + mttgl[41] + 
   mttgl[11];
   mttgl[45]=1./2.*mttgl[4];
   mttgl[41]=mttgl[45]*mttgl[41];
   mttgl[45]=mttgl[29] - mttgl[21] - 9*mttgl[20] + mttgl[31] + 
   mttgl[27];
   mttgl[43]=mttgl[43]*mttgl[45];
   mttgl[45]= - 13./2. - mttgl[14];
   mttgl[45]=3*mttgl[45] + mttgl[9];
   mttgl[36]=mttgl[45]*mttgl[36];
   mttgl[45]=mttgl[25] - 33./8. - mttgl[23];
   mttgl[36]= - 21./4.*mttgl[16] + 1./2.*mttgl[34] + 5./4.*mttgl[13] + 
   mttgl[36] + 3./2.*mttgl[33] + 3./4.*mttgl[22] - 23./4.*mttgl[19] + 
   mttgl[43] + 9./4.*mttgl[14] + 3*mttgl[45] + mttgl[41];
   mttgl[41]=1./4.*MMH;
   mttgl[36]=mttgl[41]*mttgl[36];
   mttgl[43]= - mttgl[6] + 5./2.*mttgl[12];
   mttgl[38]= - 5./2.*mttgl[7] - 23*mttgl[42] - 3./2.*mttgl[5] + 3*
   mttgl[38] + 29./2.*mttgl[17] + 13*mttgl[43] + 37*mttgl[18];
   mttgl[42]=33./2. + mttgl[14];
   mttgl[43]=mttgl[4]*mttgl[12];
   mttgl[42]= - 1./8.*mttgl[43] + 3./16.*mttgl[42] + mttgl[9];
   mttgl[42]=mttgl[11]*mttgl[42];
   mttgl[43]=mttgl[4]*mttgl[44];
   mttgl[36]=mttgl[36] + mttgl[42] - 3./4.*mttgl[24] - 3./8.*mttgl[43]
    + 1./8.*mttgl[38] + mttgl[39];
   mttgl[36]=mttgl[41]*mttgl[40]*mttgl[36];

      return mttgl[36] + mttgl[37];
}
