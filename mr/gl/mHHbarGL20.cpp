#include <HH.hpp>
namespace mr
{
  long double HH<OS>::xgl20(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<long double> armHHbarGL[56], mHHbarGLret;

    armHHbarGL[1]=double(boson);
    armHHbarGL[2]=pow(SW,-1);
    armHHbarGL[3]=pow(MMH,-1);
    armHHbarGL[4]=pow(MMW,-1);
    armHHbarGL[5]=Tsil::I2(MMH,MMH,MMH,mu2);
    armHHbarGL[6]=Tsil::I2(MMH,MMt,MMt,mu2);
    armHHbarGL[7]=Tsil::I2(0,0,MMH,mu2);
    armHHbarGL[8]=Tsil::I2(0,0,MMt,mu2);
    armHHbarGL[9]=Tsil::B(MMH,MMH,MMH,mu2);
    armHHbarGL[10]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbarGL[11]=Tsil::A(MMH,mu2);
    armHHbarGL[12]=Tsil::A(MMt,mu2);
    armHHbarGL[13]=std::real(Tsil::B(0,0,MMH,mu2));
    armHHbarGL[14]=Tsil::B(MMH,MMt,MMt,mu2);
    armHHbarGL[15]=std::real(Tsil::B(0,0,MMt,mu2));
    armHHbarGL[16]=Tsil::B(0,MMt,MMH,mu2);
    armHHbarGL[17]=Tsil::B(0,0,MMH,mu2);
    armHHbarGL[18]=Tsil::Beps(MMH,MMH,MMH,mu2);
    armHHbarGL[19]=Tsil::Beps(MMt,MMt,MMH,mu2);
    armHHbarGL[20]=Tsil::Aeps(MMH,mu2);
    armHHbarGL[21]=Tsil::Aeps(MMt,mu2);
    armHHbarGL[22]=prottttt0->Suxv(0);
    armHHbarGL[23]=protHHHHH->M(0);
    armHHbarGL[24]=protHtHtt->M(0);
    armHHbarGL[25]=protttttH->M(0);
    armHHbarGL[26]=protHHHHH->Uzxyv(0);
    armHHbarGL[27]=protHtHtt->Uzxyv(0);
    armHHbarGL[28]=protHtHtt->Uuyxv(0);
    armHHbarGL[29]=protHtHtt->Svyz(0);
    armHHbarGL[30]=prot0H0H0->M(0);
    armHHbarGL[31]=prot0t0tt->M(0);
    armHHbarGL[32]=prot0t0t0->M(0);
    armHHbarGL[33]=prot0000H->M(0);
    armHHbarGL[34]=prot0H0H0->Uuyxv(0);
    armHHbarGL[35]=prot0t0t0->Uuyxv(0);
    armHHbarGL[36]=prot0H0H0->Tyzv(0);
    armHHbarGL[37]=prot0t0t0->Tyzv(0);
    armHHbarGL[38]=1/(4*MMt - MMH);
    armHHbarGL[39]=armHHbarGL[28] + armHHbarGL[16];
    armHHbarGL[40]=armHHbarGL[18] - armHHbarGL[27];
    armHHbarGL[41]=27*armHHbarGL[9] + 21 - armHHbarGL[13];
    armHHbarGL[41]=1./4.*armHHbarGL[41] + 3*armHHbarGL[19];
    armHHbarGL[42]=1./4.*armHHbarGL[13] + armHHbarGL[9];
    armHHbarGL[42]=3*armHHbarGL[42] - 1./16.*armHHbarGL[10];
    armHHbarGL[42]=armHHbarGL[10]*armHHbarGL[42];
    armHHbarGL[39]=3./8.*armHHbarGL[17] + 1./4.*armHHbarGL[41] + 
      armHHbarGL[42] + 3*armHHbarGL[40] - 3./4.*armHHbarGL[39];
    armHHbarGL[39]=MMH*armHHbarGL[39];
    armHHbarGL[41]=5*armHHbarGL[12];
    armHHbarGL[42]=1./2.*MMH;
    armHHbarGL[43]= - armHHbarGL[10]*armHHbarGL[42];
    armHHbarGL[43]=armHHbarGL[41] + armHHbarGL[43];
    armHHbarGL[44]=1./2.*armHHbarGL[14];
    armHHbarGL[43]=armHHbarGL[43]*armHHbarGL[44];
    armHHbarGL[45]=5*armHHbarGL[10];
    armHHbarGL[46]=armHHbarGL[3]*armHHbarGL[12];
    armHHbarGL[47]= - 11*armHHbarGL[46] + 3./2. + armHHbarGL[45];
    armHHbarGL[47]=armHHbarGL[11]*armHHbarGL[47];
    armHHbarGL[48]=armHHbarGL[6] - armHHbarGL[20];
    armHHbarGL[49]=pow(armHHbarGL[12],2);
    armHHbarGL[50]=armHHbarGL[49]*armHHbarGL[3];
    armHHbarGL[51]=3./4.*armHHbarGL[10] - 9./2.*armHHbarGL[9] + 25./4.
      - armHHbarGL[13];
    armHHbarGL[51]=armHHbarGL[12]*armHHbarGL[51];
    armHHbarGL[52]=armHHbarGL[37]*armHHbarGL[42];
    armHHbarGL[39]=1./4.*armHHbarGL[47] + 13./4.*armHHbarGL[50] + 
      armHHbarGL[52] + armHHbarGL[43] + armHHbarGL[51] + 5*armHHbarGL[21]
      - armHHbarGL[22] + armHHbarGL[39] - 3./2.*armHHbarGL[48];
    armHHbarGL[43]=pow(armHHbarGL[2],4);
    armHHbarGL[39]=armHHbarGL[39]*armHHbarGL[43];
    armHHbarGL[47]=armHHbarGL[19] + 3./8. - armHHbarGL[35];
    armHHbarGL[48]=9*armHHbarGL[9];
    armHHbarGL[47]= - armHHbarGL[48] + 5./2.*armHHbarGL[47];
    armHHbarGL[51]=1 + 7*armHHbarGL[10];
    armHHbarGL[44]=armHHbarGL[51]*armHHbarGL[44];
    armHHbarGL[51]= - 1./8.*armHHbarGL[10] - armHHbarGL[48] - 5./4. - 
      armHHbarGL[13];
    armHHbarGL[51]=armHHbarGL[10]*armHHbarGL[51];
    armHHbarGL[52]=1./4.*MMH;
    armHHbarGL[53]=armHHbarGL[25]*armHHbarGL[52];
    armHHbarGL[54]= - armHHbarGL[31]*armHHbarGL[42];
    armHHbarGL[40]=armHHbarGL[54] + armHHbarGL[53] + armHHbarGL[44] + 9./
      8.*armHHbarGL[16] + 1./2.*armHHbarGL[47] + armHHbarGL[51] - 9./2.*
      armHHbarGL[40];
    armHHbarGL[44]=105./4. + armHHbarGL[10];
    armHHbarGL[44]=armHHbarGL[12]*armHHbarGL[44];
    armHHbarGL[44]=armHHbarGL[44] + 7./2.*armHHbarGL[8] - 5./2.*
      armHHbarGL[20] - 9*armHHbarGL[21] - 1./2.*armHHbarGL[22];
    armHHbarGL[41]= - armHHbarGL[14]*armHHbarGL[41];
    armHHbarGL[47]= - 4*armHHbarGL[46] - 15./8. - 2*armHHbarGL[10];
    armHHbarGL[47]=armHHbarGL[11]*armHHbarGL[47];
    armHHbarGL[41]=armHHbarGL[47] + 11./4.*armHHbarGL[6] - 2*
      armHHbarGL[50] - 17./8.*armHHbarGL[29] + 1./4.*armHHbarGL[44] + 
      armHHbarGL[41];
    armHHbarGL[41]=armHHbarGL[3]*armHHbarGL[41];
    armHHbarGL[40]=1./2.*armHHbarGL[40] + armHHbarGL[41];
    armHHbarGL[40]=armHHbarGL[40]*armHHbarGL[43];
    armHHbarGL[41]=1./4.*armHHbarGL[10];
    armHHbarGL[44]= - 1 - 7./2.*armHHbarGL[10];
    armHHbarGL[44]=armHHbarGL[44]*armHHbarGL[41];
    armHHbarGL[45]=armHHbarGL[45] + 1;
    armHHbarGL[47]= - armHHbarGL[14]*armHHbarGL[45];
    armHHbarGL[50]=2*armHHbarGL[21] + armHHbarGL[20];
    armHHbarGL[50]=9./2.*armHHbarGL[12] + 2*armHHbarGL[50] - 1./2.*
      armHHbarGL[8];
    armHHbarGL[50]=armHHbarGL[3]*armHHbarGL[50];
    armHHbarGL[44]=armHHbarGL[50] - 3./8.*armHHbarGL[37] + 
      armHHbarGL[47] + 2*armHHbarGL[28] + 5./16.*armHHbarGL[16] + 
      armHHbarGL[44] - 13./4.*armHHbarGL[19] + 6 + 5./4.*armHHbarGL[35];
    armHHbarGL[44]=armHHbarGL[3]*armHHbarGL[44];
    armHHbarGL[47]= - 2 - armHHbarGL[10];
    armHHbarGL[47]=armHHbarGL[11]*armHHbarGL[47];
    armHHbarGL[47]=armHHbarGL[47] - armHHbarGL[6];
    armHHbarGL[47]=2*armHHbarGL[47];
    armHHbarGL[47]=armHHbarGL[47]*pow(armHHbarGL[3],2);
    armHHbarGL[44]=armHHbarGL[44] - 1./2.*armHHbarGL[32] + 
      armHHbarGL[25] + armHHbarGL[47];
    armHHbarGL[44]=armHHbarGL[44]*armHHbarGL[43];
    armHHbarGL[47]=armHHbarGL[43]*MMt;
    armHHbarGL[50]=armHHbarGL[47]*armHHbarGL[3];
    armHHbarGL[51]= - armHHbarGL[37] - 9 - armHHbarGL[16];
    armHHbarGL[51]=armHHbarGL[3]*armHHbarGL[51];
    armHHbarGL[51]= - 2*armHHbarGL[25] + 1./2.*armHHbarGL[51];
    armHHbarGL[51]=armHHbarGL[51]*armHHbarGL[50];
    armHHbarGL[44]=armHHbarGL[44] + armHHbarGL[51];
    armHHbarGL[44]=MMt*armHHbarGL[44];
    armHHbarGL[40]=armHHbarGL[40] + armHHbarGL[44];
    armHHbarGL[40]=MMt*armHHbarGL[40];
    armHHbarGL[39]=1./2.*armHHbarGL[39] + armHHbarGL[40];
    armHHbarGL[39]=MMt*armHHbarGL[39];
    armHHbarGL[40]=1./2.*armHHbarGL[10];
    armHHbarGL[44]=armHHbarGL[40] - 1;
    armHHbarGL[41]=armHHbarGL[44]*armHHbarGL[41];
    armHHbarGL[51]= - 7 + armHHbarGL[13];
    armHHbarGL[51]=armHHbarGL[13]*armHHbarGL[51];
    armHHbarGL[51]= - 25 + armHHbarGL[51];
    armHHbarGL[53]=3*armHHbarGL[13];
    armHHbarGL[54]=21./4.*armHHbarGL[9] + 13./4. + armHHbarGL[53];
    armHHbarGL[54]=armHHbarGL[9]*armHHbarGL[54];
    armHHbarGL[55]=armHHbarGL[30] + 9./2.*armHHbarGL[23];
    armHHbarGL[55]=3*armHHbarGL[55] + 1./2.*armHHbarGL[33];
    armHHbarGL[55]=MMH*armHHbarGL[55];
    armHHbarGL[41]=armHHbarGL[55] + 3./2.*armHHbarGL[17] + 
      armHHbarGL[41] - 9./2.*armHHbarGL[26] + 1./4.*armHHbarGL[51] + 3*
      armHHbarGL[54];
    armHHbarGL[41]=armHHbarGL[41]*armHHbarGL[52];
    armHHbarGL[48]= - armHHbarGL[10] - armHHbarGL[48] - 5 + 
      armHHbarGL[53];
    armHHbarGL[48]=armHHbarGL[12]*armHHbarGL[48];
    armHHbarGL[41]=armHHbarGL[41] + 1./4.*armHHbarGL[48] + 1./8.*
      armHHbarGL[20] + 3./8.*armHHbarGL[5] + armHHbarGL[21];
    armHHbarGL[41]=MMH*armHHbarGL[41];
    armHHbarGL[48]=pow(Pi,2);
    armHHbarGL[48]=11./8.*armHHbarGL[36] + 9./4.*armHHbarGL[18] + 5./16.
      *armHHbarGL[48];
    armHHbarGL[51]=pow(MMH,2);
    armHHbarGL[48]=armHHbarGL[51]*armHHbarGL[48];
    armHHbarGL[44]=armHHbarGL[10]*armHHbarGL[44];
    armHHbarGL[44]=1./2. + armHHbarGL[44];
    armHHbarGL[44]=armHHbarGL[44]*armHHbarGL[52];
    armHHbarGL[52]= - 1 + armHHbarGL[10];
    armHHbarGL[52]=armHHbarGL[12]*armHHbarGL[52];
    armHHbarGL[44]=armHHbarGL[52] + armHHbarGL[44];
    armHHbarGL[44]=armHHbarGL[44]*armHHbarGL[42];
    armHHbarGL[44]=armHHbarGL[49] + armHHbarGL[44];
    armHHbarGL[42]=armHHbarGL[38]*armHHbarGL[44]*armHHbarGL[42];
    armHHbarGL[41]=armHHbarGL[42] - 1./2.*armHHbarGL[49] + 
      armHHbarGL[41] + armHHbarGL[48];
    armHHbarGL[42]=45./16.*armHHbarGL[9] - 1 + 13./16.*armHHbarGL[13];
    armHHbarGL[42]=MMH*armHHbarGL[42];
    armHHbarGL[42]=51./32.*armHHbarGL[11] - 3./4.*armHHbarGL[12] + 
      armHHbarGL[42];
    armHHbarGL[42]=armHHbarGL[11]*armHHbarGL[42];
    armHHbarGL[44]=armHHbarGL[34]*armHHbarGL[51];
    armHHbarGL[41]= - 9./16.*armHHbarGL[44] + 1./2.*armHHbarGL[41] + 
      armHHbarGL[42];
    armHHbarGL[41]=armHHbarGL[41]*armHHbarGL[43];
    armHHbarGL[42]=armHHbarGL[43]*MMH;
    armHHbarGL[44]=armHHbarGL[7]*armHHbarGL[42];
    armHHbarGL[41]=armHHbarGL[41] + 1./16.*armHHbarGL[44];
    armHHbarGL[39]=1./2.*armHHbarGL[41] + armHHbarGL[39];
    armHHbarGL[39]=armHHbarGL[39]*pow(armHHbarGL[4],2);
    armHHbarGL[40]=armHHbarGL[40] - 5*armHHbarGL[46];
    armHHbarGL[40]=armHHbarGL[40]*armHHbarGL[43];
    armHHbarGL[41]= - armHHbarGL[45]*armHHbarGL[50];
    armHHbarGL[40]=armHHbarGL[40] + armHHbarGL[41];
    armHHbarGL[40]=armHHbarGL[15]*armHHbarGL[40];
    armHHbarGL[41]=3./4.*armHHbarGL[42] - 2*armHHbarGL[47];
    armHHbarGL[41]=armHHbarGL[24]*armHHbarGL[41];
    armHHbarGL[40]=3*armHHbarGL[41] + 1./4.*armHHbarGL[40];
    armHHbarGL[41]=armHHbarGL[4]*MMt;
    armHHbarGL[40]=armHHbarGL[40]*pow(armHHbarGL[41],2);
    armHHbarGL[39]=armHHbarGL[39] + armHHbarGL[40];

    mHHbarGLret = 3*armHHbarGL[39]*armHHbarGL[1];
    return mHHbarGLret.real();
  }
} // namespace mr
