#include <tt.hpp>
namespace mr
{
  double tt<OS>::ygl20(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryuttGL[49], yuttGLret;

    aryuttGL[1]=double(boson);
    aryuttGL[2]=pow(SW,-1);
    aryuttGL[3]=pow(MMH,-1);
    aryuttGL[4]=pow(MMW,-1);
    aryuttGL[5]=pow(MMt,-1);
    aryuttGL[6]=Tsil::I2(MMH,MMH,MMH,mu2);
    aryuttGL[7]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryuttGL[8]=Tsil::I2(0,MMH,MMt,mu2);
    aryuttGL[9]=Tsil::I2(0,0,MMH,mu2);
    aryuttGL[10]=Tsil::I2(0,0,MMt,mu2);
    aryuttGL[11]=Tsil::B(MMH,MMH,MMH,mu2);
    aryuttGL[12]=Tsil::B(MMH,MMt,MMt,mu2);
    aryuttGL[13]=Tsil::A(MMH,mu2);
    aryuttGL[14]=Tsil::A(MMt,mu2);
    aryuttGL[15]=Tsil::B(MMt,MMt,MMH,mu2);
    aryuttGL[16]=std::real(Tsil::B(0,0,MMH,mu2));
    aryuttGL[17]=std::real(Tsil::B(0,0,MMt,mu2));
    aryuttGL[18]=Tsil::B(0,MMH,MMt,mu2);
    aryuttGL[19]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryuttGL[20]=Tsil::Aeps(MMH,mu2);
    aryuttGL[21]=Tsil::Aeps(MMt,mu2);
    aryuttGL[22]=prot0ttHt->Tuxv(0);
    aryuttGL[23]=protHHttH->M(0);
    aryuttGL[24]=protHttHt->M(0);
    aryuttGL[25]=protHHttH->Uxzuv(0);
    aryuttGL[26]=protHHttH->Suxv(0);
    aryuttGL[27]=protHHttH->Uuyxv(0);
    aryuttGL[28]=protWt000->Tyzv(0);
    aryuttGL[29]=protH0tt0->M(0);
    aryuttGL[30]=protH0t00->M(0);
    aryuttGL[31]=prot0Htt0->M(0);
    aryuttGL[32]=prot0H0t0->M(0);
    aryuttGL[33]=prot00ttH->M(0);
    aryuttGL[34]=protH0t00->Uxzuv(0);
    aryuttGL[35]=protH0t00->Uzxyv(0);
    aryuttGL[36]=protH0tt0->Uuyxv(0);
    aryuttGL[37]=protH0t00->Txuv(0);
    aryuttGL[38]=3*aryuttGL[15];
    aryuttGL[39]=1 - aryuttGL[12];
    aryuttGL[39]=aryuttGL[38]*aryuttGL[39];
    aryuttGL[39]=aryuttGL[39] - 111./16.;
    aryuttGL[39]=aryuttGL[3]*aryuttGL[39];
    aryuttGL[40]=3*aryuttGL[3];
    aryuttGL[41]= - aryuttGL[12]*aryuttGL[40];
    aryuttGL[39]=aryuttGL[41] + aryuttGL[24] + aryuttGL[39];
    aryuttGL[39]=MMt*aryuttGL[39];
    aryuttGL[41]=aryuttGL[7]*aryuttGL[3];
    aryuttGL[41]= - 1073./64. - 9*aryuttGL[41];
    aryuttGL[40]=aryuttGL[10]*aryuttGL[40];
    aryuttGL[42]=aryuttGL[8]*aryuttGL[3];
    aryuttGL[43]=aryuttGL[13]*pow(aryuttGL[3],2);
    aryuttGL[43]= - aryuttGL[3] - 21./2.*aryuttGL[43];
    aryuttGL[43]=aryuttGL[13]*aryuttGL[43];
    aryuttGL[39]=aryuttGL[39] - 115./16.*aryuttGL[22] - 1./2.*
      aryuttGL[25] + 3./8.*aryuttGL[43] - 9./8.*aryuttGL[42] + 1./4.*
      aryuttGL[41] + aryuttGL[40];
    aryuttGL[40]=aryuttGL[14]*aryuttGL[3];
    aryuttGL[41]=aryuttGL[13]*aryuttGL[3];
    aryuttGL[42]=3./2.*aryuttGL[12] - 3*aryuttGL[40] + 11./16.*
      aryuttGL[17] - 45./4. + aryuttGL[41];
    aryuttGL[43]=1./2.*aryuttGL[12];
    aryuttGL[42]=aryuttGL[42]*aryuttGL[43];
    aryuttGL[44]=5./4.*aryuttGL[12] - aryuttGL[40] - 1./2. + 
      aryuttGL[41];
    aryuttGL[45]=3./2.*aryuttGL[15];
    aryuttGL[44]=aryuttGL[44]*aryuttGL[45];
    aryuttGL[46]= - 185./32.*aryuttGL[20] - 6*aryuttGL[21];
    aryuttGL[46]=aryuttGL[3]*aryuttGL[46];
    aryuttGL[47]=3./8.*aryuttGL[17] - 39./8. + aryuttGL[41];
    aryuttGL[47]=aryuttGL[17]*aryuttGL[47];
    aryuttGL[48]=aryuttGL[30] + aryuttGL[32];
    aryuttGL[48]=3*aryuttGL[23] - aryuttGL[24] + 1./4.*aryuttGL[48];
    aryuttGL[48]=MMH*aryuttGL[48];
    aryuttGL[39]=1./4.*aryuttGL[48] + aryuttGL[44] - 3./128.*
      aryuttGL[18] + 75./64.*aryuttGL[28] + aryuttGL[42] + 21./16.*
      aryuttGL[40] + 13./32.*aryuttGL[19] + 1./16.*aryuttGL[47] - 5./32.*
      aryuttGL[34] + aryuttGL[46] + 1./2.*aryuttGL[39];
    aryuttGL[39]=MMt*aryuttGL[39];
    aryuttGL[42]=1./4. - aryuttGL[12];
    aryuttGL[42]=aryuttGL[42]*aryuttGL[45];
    aryuttGL[44]=3*aryuttGL[16] + 9*aryuttGL[11];
    aryuttGL[45]=361./4. - aryuttGL[44];
    aryuttGL[45]= - 9./2.*aryuttGL[27] + 35./4.*aryuttGL[19] + 1./8.*
      aryuttGL[17] + 5./4.*aryuttGL[34] + 1./2.*aryuttGL[45] + 35*
      aryuttGL[22];
    aryuttGL[46]=13 + aryuttGL[11];
    aryuttGL[46]=3*aryuttGL[46] + aryuttGL[16];
    aryuttGL[46]=3./8.*aryuttGL[46] - aryuttGL[12];
    aryuttGL[46]=aryuttGL[12]*aryuttGL[46];
    aryuttGL[47]= - 9*aryuttGL[23] + aryuttGL[33] + aryuttGL[31] + 
      aryuttGL[29] - aryuttGL[24];
    aryuttGL[47]=MMH*aryuttGL[47];
    aryuttGL[42]=1./8.*aryuttGL[47] + aryuttGL[42] - 3./32.*aryuttGL[18]
      + 5./8.*aryuttGL[28] + 1./4.*aryuttGL[45] + aryuttGL[46];
    aryuttGL[42]=MMH*aryuttGL[42];
    aryuttGL[45]=1./4.*aryuttGL[13];
    aryuttGL[46]= - 261./2. + 61*aryuttGL[41];
    aryuttGL[46]=aryuttGL[46]*aryuttGL[45];
    aryuttGL[41]= - 3*aryuttGL[17] - 399./4. + 125*aryuttGL[41];
    aryuttGL[40]=1./2.*aryuttGL[41] - 93*aryuttGL[40];
    aryuttGL[40]=aryuttGL[14]*aryuttGL[40];
    aryuttGL[40]=aryuttGL[40] + 85./4.*aryuttGL[20] + aryuttGL[46] + 27./
      2.*aryuttGL[8] + 5./2.*aryuttGL[10] + 17./4.*aryuttGL[26] + 15*
      aryuttGL[7];
    aryuttGL[41]= - aryuttGL[13] + aryuttGL[14];
    aryuttGL[38]=aryuttGL[41]*aryuttGL[38];
    aryuttGL[41]= - 1./2.*aryuttGL[13] + 3*aryuttGL[14];
    aryuttGL[41]=aryuttGL[12]*aryuttGL[41];
    aryuttGL[38]=aryuttGL[38] + 1./2.*aryuttGL[40] + 3*aryuttGL[41];
    aryuttGL[38]=1./2.*aryuttGL[38] + aryuttGL[42];
    aryuttGL[38]=1./4.*aryuttGL[38] + aryuttGL[39];
    aryuttGL[38]=MMt*aryuttGL[38];
    aryuttGL[39]= - 1 - aryuttGL[18];
    aryuttGL[39]=MMH*aryuttGL[39];
    aryuttGL[39]=aryuttGL[39] + aryuttGL[14] - aryuttGL[9] + 
      aryuttGL[13] + aryuttGL[21];
    aryuttGL[40]=1./2.*aryuttGL[5];
    aryuttGL[39]=aryuttGL[40]*aryuttGL[39];
    aryuttGL[41]=aryuttGL[12] - 61./4. - aryuttGL[44];
    aryuttGL[41]=aryuttGL[41]*aryuttGL[43];
    aryuttGL[42]= - 21./2.*aryuttGL[19] - 23./2.*aryuttGL[22] + 3./2.*
      aryuttGL[25] + 3./4.*aryuttGL[16] + 9./4.*aryuttGL[11] - 571./16. + 
      aryuttGL[36];
    aryuttGL[39]=5./4.*aryuttGL[18] + aryuttGL[41] + 243./4.*S2 + 1./2.*
      aryuttGL[42] + 3*aryuttGL[27] + aryuttGL[39];
    aryuttGL[39]=MMH*aryuttGL[39];
    aryuttGL[41]=3 + aryuttGL[11];
    aryuttGL[41]=3*aryuttGL[41] + aryuttGL[16];
    aryuttGL[42]= - aryuttGL[13] + 3./2.*aryuttGL[14];
    aryuttGL[42]=aryuttGL[5]*aryuttGL[42];
    aryuttGL[41]=1./4.*aryuttGL[42] + 3./4.*aryuttGL[41];
    aryuttGL[41]=aryuttGL[14]*aryuttGL[41];
    aryuttGL[42]=131./4. - aryuttGL[44];
    aryuttGL[42]=aryuttGL[42]*aryuttGL[45];
    aryuttGL[44]=aryuttGL[45] - 7*aryuttGL[14];
    aryuttGL[43]=aryuttGL[44]*aryuttGL[43];
    aryuttGL[39]=1./2.*aryuttGL[39] + aryuttGL[43] - 3./8.*aryuttGL[6]
      - 5./8.*aryuttGL[9] + 17./8.*aryuttGL[20] + 31./4.*aryuttGL[21] + 
      aryuttGL[42] - aryuttGL[7] - 9./4.*aryuttGL[8] + aryuttGL[41];
    aryuttGL[39]=MMH*aryuttGL[39];
    aryuttGL[41]=pow(aryuttGL[13],2);
    aryuttGL[42]= - 9*aryuttGL[13] + 25*aryuttGL[14];
    aryuttGL[42]=aryuttGL[14]*aryuttGL[42];
    aryuttGL[41]= - 19./4.*aryuttGL[41] + aryuttGL[42];
    aryuttGL[39]=1./2.*aryuttGL[41] + aryuttGL[39];
    aryuttGL[40]= - pow(MMH,3)*aryuttGL[40];
    aryuttGL[41]=3./4.*MMH + MMt;
    aryuttGL[41]=MMt*aryuttGL[41];
    aryuttGL[40]=aryuttGL[40] + aryuttGL[41];
    aryuttGL[40]=aryuttGL[37]*aryuttGL[40];
    aryuttGL[41]=pow(MMH,2);
    aryuttGL[42]=MMt*MMH;
    aryuttGL[42]=3*aryuttGL[41] - 11*aryuttGL[42];
    aryuttGL[42]=aryuttGL[35]*aryuttGL[42];
    aryuttGL[43]=MMt*aryuttGL[3];
    aryuttGL[43]= - 1 - 29./16.*aryuttGL[43];
    aryuttGL[43]=MMt*aryuttGL[43];
    aryuttGL[43]=9./16.*MMH + aryuttGL[43];
    aryuttGL[43]=MMt*aryuttGL[43];
    aryuttGL[41]= - 1./12.*aryuttGL[41] + aryuttGL[43];
    aryuttGL[41]=aryuttGL[41]*pow(Pi,2);
    aryuttGL[38]=1./4.*aryuttGL[41] + 1./32.*aryuttGL[42] + 1./16.*
      aryuttGL[40] + 1./8.*aryuttGL[39] + aryuttGL[38];

    yuttGLret = aryuttGL[38]*pow(aryuttGL[4],2)*pow(aryuttGL[2],4)*
      aryuttGL[1];
    return yuttGLret.real();
  }
} // namespace mr
