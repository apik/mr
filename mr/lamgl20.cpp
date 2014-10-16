#include <HH.hpp>
std::complex<long double> HH::mygl20(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mlamgl[50];

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
    mlamgl[27]=prot0H0H0->M(0);
    mlamgl[28]=prot0t0tt->M(0);
    mlamgl[29]=prot0t0t0->M(0);
    mlamgl[30]=prot0000H->M(0);
    mlamgl[31]=prot0H0H0->Uuyxv(0);
    mlamgl[32]=prot0t0t0->Uuyxv(0);
    mlamgl[33]=prot0H0H0->Tyzv(0);
    mlamgl[34]=prot0t0t0->Tyzv(0);
    mlamgl[35]=1/(4*MMt - MMH);
   mlamgl[36]=1./4.*mlamgl[8];
   mlamgl[37]=55 - mlamgl[8];
   mlamgl[37]=mlamgl[37]*mlamgl[36];
   mlamgl[38]=1./2.*mlamgl[13];
   mlamgl[39]=mlamgl[8] - 1./2.;
   mlamgl[40]= - mlamgl[9] + mlamgl[39];
   mlamgl[40]=mlamgl[40]*mlamgl[38];
   mlamgl[41]=3*mlamgl[8];
   mlamgl[42]=1 - mlamgl[41];
   mlamgl[42]=mlamgl[11]*mlamgl[42];
   mlamgl[42]= - mlamgl[32] + mlamgl[42] + mlamgl[15];
   mlamgl[43]=9*mlamgl[14];
   mlamgl[44]=7*mlamgl[8];
   mlamgl[45]= - 3 + mlamgl[44];
   mlamgl[45]=mlamgl[10]*mlamgl[45];
   mlamgl[41]=37 - mlamgl[41];
   mlamgl[41]=1./2.*mlamgl[41] - mlamgl[10];
   mlamgl[41]=mlamgl[9]*mlamgl[41];
   mlamgl[37]=mlamgl[40] - mlamgl[43] + mlamgl[41] + mlamgl[45] - 23 + 
   mlamgl[37] + 5./2.*mlamgl[42];
   mlamgl[40]=mlamgl[16]*mlamgl[11];
   mlamgl[41]=3*mlamgl[40];
   mlamgl[42]=pow(Pi,2);
   mlamgl[44]= - mlamgl[44] - 3*mlamgl[9];
   mlamgl[44]=mlamgl[7]*mlamgl[44];
   mlamgl[45]=1./2.*mlamgl[9];
   mlamgl[46]= - 9*mlamgl[11] + mlamgl[45];
   mlamgl[46]=mlamgl[17]*mlamgl[46];
   mlamgl[37]=3./2.*mlamgl[46] + 27*mlamgl[23] - mlamgl[41] + 27./4.*
   mlamgl[12] + 27./2.*mlamgl[44] + 3*mlamgl[37] - 17./16.*mlamgl[42];
   mlamgl[44]= - 1./4.*mlamgl[18] + mlamgl[5] - 5./4.*mlamgl[26];
   mlamgl[44]= - 9./4.*mlamgl[6] + 3*mlamgl[44];
   mlamgl[44]=mlamgl[2]*mlamgl[44];
   mlamgl[46]=9*mlamgl[20] + 1./2.*mlamgl[21] - mlamgl[28];
   mlamgl[46]=MMH*mlamgl[46];
   mlamgl[37]=3./2.*mlamgl[46] - 9*mlamgl[25] + 1./2.*mlamgl[37] + 
   mlamgl[44];
   mlamgl[44]= - mlamgl[45] + 1 - 5./2.*mlamgl[8];
   mlamgl[38]=mlamgl[44]*mlamgl[38];
   mlamgl[44]=mlamgl[17]*mlamgl[9];
   mlamgl[46]=pow(mlamgl[8],2);
   mlamgl[46]= - 7*mlamgl[46] + 23./2. - 3*mlamgl[34];
   mlamgl[46]=1./2.*mlamgl[46] + 5*mlamgl[32];
   mlamgl[47]=2 - 5*mlamgl[8];
   mlamgl[47]=mlamgl[10]*mlamgl[47];
   mlamgl[48]= - 2 - mlamgl[8];
   mlamgl[48]=mlamgl[11]*mlamgl[48];
   mlamgl[49]=29./4. + 13*mlamgl[8];
   mlamgl[49]=1./4.*mlamgl[49] - mlamgl[10];
   mlamgl[49]=mlamgl[9]*mlamgl[49];
   mlamgl[38]=mlamgl[44] + mlamgl[38] + mlamgl[49] - 13./4.*mlamgl[15]
    + 2*mlamgl[48] + 1./4.*mlamgl[46] + mlamgl[47];
   mlamgl[46]= - mlamgl[5] + 2*mlamgl[26];
   mlamgl[46]=mlamgl[2]*mlamgl[46];
   mlamgl[46]=mlamgl[46] + mlamgl[24];
   mlamgl[38]= - 15./2.*mlamgl[25] + mlamgl[41] + 15./16.*mlamgl[12] + 
   7./8.*mlamgl[42] + 3*mlamgl[38] + 6*mlamgl[46];
   mlamgl[38]=mlamgl[2]*mlamgl[38];
   mlamgl[41]= - mlamgl[12] + mlamgl[44] + 51./2. - mlamgl[34];
   mlamgl[41]=1./8.*mlamgl[42] - 9*mlamgl[9] + 1./2.*mlamgl[41];
   mlamgl[41]=mlamgl[2]*mlamgl[41];
   mlamgl[41]= - 2*mlamgl[21] + mlamgl[41];
   mlamgl[41]=mlamgl[2]*mlamgl[41];
   mlamgl[46]=mlamgl[25]*pow(mlamgl[2],2);
   mlamgl[41]=mlamgl[41] + 4*mlamgl[46];
   mlamgl[41]=MMt*mlamgl[41];
   mlamgl[41]=mlamgl[41] - 6*mlamgl[20] + mlamgl[21] - 1./2.*mlamgl[29]
   ;
   mlamgl[38]=mlamgl[38] + 3*mlamgl[41];
   mlamgl[38]=MMt*mlamgl[38];
   mlamgl[37]=1./2.*mlamgl[37] + mlamgl[38];
   mlamgl[38]=pow(mlamgl[1],4)*pow(mlamgl[3],2);
   mlamgl[37]=MMt*mlamgl[38]*mlamgl[37];
   mlamgl[41]=1./2.*mlamgl[8];
   mlamgl[46]=11 - mlamgl[41];
   mlamgl[46]=mlamgl[46]*mlamgl[36];
   mlamgl[47]=1./2.*mlamgl[10];
   mlamgl[39]= - mlamgl[39]*mlamgl[47];
   mlamgl[39]=mlamgl[39] + mlamgl[46] + 43./8. + mlamgl[34];
   mlamgl[46]=mlamgl[47] - 3 - mlamgl[36];
   mlamgl[45]=mlamgl[46]*mlamgl[45];
   mlamgl[46]= - 1./4. - mlamgl[8];
   mlamgl[46]=mlamgl[11]*mlamgl[46];
   mlamgl[39]=3*mlamgl[14] + mlamgl[45] + 3./4.*mlamgl[15] + 1./2.*
   mlamgl[39] + mlamgl[46];
   mlamgl[45]=9*mlamgl[7];
   mlamgl[46]= - 3./4.*mlamgl[9] + 3./2. + mlamgl[8];
   mlamgl[46]=mlamgl[46]*mlamgl[45];
   mlamgl[47]=mlamgl[24] + mlamgl[12];
   mlamgl[48]=mlamgl[17]*mlamgl[11];
   mlamgl[39]=15./8.*mlamgl[48] - 9*mlamgl[23] + 9./16.*mlamgl[40] + 
   mlamgl[46] + 3*mlamgl[39] + 1./16.*mlamgl[42] - 9./4.*mlamgl[47];
   mlamgl[39]=MMH*mlamgl[39];
   mlamgl[46]= - 3./4.*mlamgl[5] - mlamgl[18];
   mlamgl[39]=27./8.*mlamgl[6] + 3*mlamgl[46] + mlamgl[39];
   mlamgl[39]=mlamgl[39]*mlamgl[38];
   mlamgl[37]=1./2.*mlamgl[39] + mlamgl[37];
   mlamgl[37]=MMt*mlamgl[37];
   mlamgl[39]=mlamgl[41] - 1;
   mlamgl[39]=mlamgl[39]*mlamgl[9];
   mlamgl[36]=mlamgl[36] - 1;
   mlamgl[36]=mlamgl[36]*mlamgl[8];
   mlamgl[36]=mlamgl[39] + mlamgl[36];
   mlamgl[39]=mlamgl[31] + mlamgl[22];
   mlamgl[39]= - 51./4.*mlamgl[11] - 131./4. + mlamgl[36] - 9*
   mlamgl[39];
   mlamgl[36]=mlamgl[36] + 1./4.*mlamgl[44] + 1;
   mlamgl[36]=mlamgl[35]*mlamgl[36];
   mlamgl[36]=27*mlamgl[19] + mlamgl[30] + mlamgl[36];
   mlamgl[36]=1./2.*mlamgl[36] + 3*mlamgl[27];
   mlamgl[36]=MMH*mlamgl[36];
   mlamgl[36]=mlamgl[36] + 1./2.*mlamgl[39] + mlamgl[43];
   mlamgl[39]=67./8. + 3*mlamgl[7];
   mlamgl[39]=mlamgl[39]*mlamgl[45];
   mlamgl[36]=33./2.*mlamgl[33] + 3./8.*mlamgl[44] + 39./2.*mlamgl[40]
    + 243./2.*S2 + mlamgl[39] - 23./12.*mlamgl[42] + 3*mlamgl[36];
   mlamgl[36]=MMH*mlamgl[36];
   mlamgl[36]=9*mlamgl[5] + mlamgl[36] + 9./2.*mlamgl[4] - 9*mlamgl[6];
   mlamgl[36]=MMH*mlamgl[38]*mlamgl[36];

      return 1./16.*mlamgl[36] + mlamgl[37];
}
