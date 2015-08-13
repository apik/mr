#include <tt.hpp>
long double tt<OS>::xgl20(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbarGL[50], mttbarGLret;

    armttbarGL[1]=double(boson);
    armttbarGL[2]=pow(SW,-1);
    armttbarGL[3]=pow(MMH,-1);
    armttbarGL[4]=pow(MMW,-1);
    armttbarGL[5]=pow(MMt,-1);
    armttbarGL[6]=Tsil::I2(MMH,MMH,MMH,mu2);
    armttbarGL[7]=Tsil::I2(MMH,MMt,MMt,mu2);
    armttbarGL[8]=Tsil::I2(0,0,MMH,mu2);
    armttbarGL[9]=Tsil::I2(0,0,MMt,mu2);
    armttbarGL[10]=Tsil::B(MMH,MMH,MMH,mu2);
    armttbarGL[11]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbarGL[12]=Tsil::A(MMH,mu2);
    armttbarGL[13]=Tsil::A(MMt,mu2);
    armttbarGL[14]=Tsil::B(MMt,MMt,MMH,mu2);
    armttbarGL[15]=std::real(Tsil::B(0,0,MMH,mu2));
    armttbarGL[16]=std::real(Tsil::B(0,0,MMt,mu2));
    armttbarGL[17]=Tsil::B(0,MMH,MMt,mu2);
    armttbarGL[18]=Tsil::Beps(MMH,MMt,MMt,mu2);
    armttbarGL[19]=Tsil::Aeps(MMH,mu2);
    armttbarGL[20]=Tsil::Aeps(MMt,mu2);
    armttbarGL[21]=prot0ttHt->Tuxv(0);
    armttbarGL[22]=protHHttH->M(0);
    armttbarGL[23]=protHttHt->M(0);
    armttbarGL[24]=protHHttH->Uxzuv(0);
    armttbarGL[25]=protHHttH->Suxv(0);
    armttbarGL[26]=protHHttH->Uuyxv(0);
    armttbarGL[27]=protWt000->Tyzv(0);
    armttbarGL[28]=protH0tt0->M(0);
    armttbarGL[29]=protH0t00->M(0);
    armttbarGL[30]=prot0Htt0->M(0);
    armttbarGL[31]=prot0H0t0->M(0);
    armttbarGL[32]=prot00ttH->M(0);
    armttbarGL[33]=protH0t00->Uxzuv(0);
    armttbarGL[34]=protH0t00->Uzxyv(0);
    armttbarGL[35]=protH0tt0->Uuyxv(0);
    armttbarGL[36]=protH0t00->Txuv(0);
   armttbarGL[37]=3*armttbarGL[11];
   armttbarGL[38]= - 19./4. + armttbarGL[11];
   armttbarGL[38]=armttbarGL[38]*armttbarGL[37];
   armttbarGL[39]=1./8.*armttbarGL[16];
   armttbarGL[40]=11*armttbarGL[11];
   armttbarGL[41]=3./4.*armttbarGL[16] + 3./4. + armttbarGL[40];
   armttbarGL[41]=armttbarGL[41]*armttbarGL[39];
   armttbarGL[42]=3./2.*armttbarGL[14];
   armttbarGL[43]= - 7 + 5*armttbarGL[11];
   armttbarGL[43]=armttbarGL[43]*armttbarGL[42];
   armttbarGL[44]=armttbarGL[31]*MMH;
   armttbarGL[44]=armttbarGL[44] + armttbarGL[36];
   armttbarGL[45]=3./32.*armttbarGL[17];
   armttbarGL[46]=1./4.*armttbarGL[29] - armttbarGL[23];
   armttbarGL[46]=MMH*armttbarGL[46];
   armttbarGL[47]=pow(Pi,2);
   armttbarGL[38]=13./8.*armttbarGL[18] + 3./32.*armttbarGL[47] + 
   armttbarGL[46] + armttbarGL[43] - armttbarGL[45] + 75./16.*
   armttbarGL[27] - 115./8.*armttbarGL[21] + armttbarGL[41] + 
   armttbarGL[38] - 473./128. - armttbarGL[24] + 1./4.*armttbarGL[44];
   armttbarGL[38]= - 5./32.*armttbarGL[33] + 1./4.*armttbarGL[38];
   armttbarGL[41]=pow(armttbarGL[2],4);
   armttbarGL[38]=armttbarGL[41]*armttbarGL[38];
   armttbarGL[43]=armttbarGL[14]*armttbarGL[13];
   armttbarGL[44]=3./2.*armttbarGL[7] + armttbarGL[9];
   armttbarGL[46]= - armttbarGL[37] - 5./8.*armttbarGL[16];
   armttbarGL[46]=armttbarGL[13]*armttbarGL[46];
   armttbarGL[44]=armttbarGL[43] + armttbarGL[46] + 1./2.*
   armttbarGL[44] - 2*armttbarGL[20];
   armttbarGL[39]= - armttbarGL[42] + armttbarGL[39] - 33./8. + 
   armttbarGL[11];
   armttbarGL[39]=armttbarGL[12]*armttbarGL[39];
   armttbarGL[39]=3*armttbarGL[44] + 1./2.*armttbarGL[39];
   armttbarGL[39]=armttbarGL[39]*armttbarGL[41];
   armttbarGL[44]=pow(armttbarGL[13],2);
   armttbarGL[46]=armttbarGL[41]*armttbarGL[3];
   armttbarGL[47]=armttbarGL[44]*armttbarGL[46];
   armttbarGL[39]=armttbarGL[39] - 9./2.*armttbarGL[47];
   armttbarGL[39]=armttbarGL[3]*armttbarGL[39];
   armttbarGL[47]=1./2.*armttbarGL[16];
   armttbarGL[48]=armttbarGL[11] - 1;
   armttbarGL[49]= - armttbarGL[14]*armttbarGL[48];
   armttbarGL[37]=armttbarGL[49] - armttbarGL[47] + 5./2. - 
   armttbarGL[37];
   armttbarGL[49]=1./2.*armttbarGL[41];
   armttbarGL[37]=armttbarGL[37]*armttbarGL[49];
   armttbarGL[43]=armttbarGL[43]*armttbarGL[46];
   armttbarGL[37]=armttbarGL[37] - 6*armttbarGL[43];
   armttbarGL[37]=armttbarGL[3]*armttbarGL[37];
   armttbarGL[43]=armttbarGL[23]*armttbarGL[49];
   armttbarGL[37]=armttbarGL[43] + 3*armttbarGL[37];
   armttbarGL[37]=armttbarGL[37]*MMt;
   armttbarGL[43]=armttbarGL[19]*armttbarGL[46];
   armttbarGL[37]=armttbarGL[37] - 185./32.*armttbarGL[43] + 
   armttbarGL[39] + armttbarGL[38];
   armttbarGL[38]=pow(armttbarGL[4],2);
   armttbarGL[37]=MMt*armttbarGL[38]*armttbarGL[37];
   armttbarGL[39]=3./8.*armttbarGL[15] + 9./8.*armttbarGL[10];
   armttbarGL[39]=armttbarGL[48]*armttbarGL[39];
   armttbarGL[43]=armttbarGL[11] - 3./2.;
   armttbarGL[42]= - armttbarGL[43]*armttbarGL[42];
   armttbarGL[46]=211./16. - armttbarGL[11];
   armttbarGL[46]=armttbarGL[11]*armttbarGL[46];
   armttbarGL[48]= - armttbarGL[23] + armttbarGL[28] + armttbarGL[32]
    + armttbarGL[30];
   armttbarGL[48]=MMH*armttbarGL[48];
   armttbarGL[39]=1./8.*armttbarGL[48] + armttbarGL[42] + 3./16.*
   armttbarGL[36] - armttbarGL[45] + 5./8.*armttbarGL[27] + 35./4.*
   armttbarGL[21] - 9./8.*armttbarGL[26] + 89./16. + armttbarGL[46] + 
   armttbarGL[39] + 5./16.*armttbarGL[33] - 11./8.*armttbarGL[34] + 35./
   16.*armttbarGL[18];
   armttbarGL[39]=MMH*armttbarGL[39];
   armttbarGL[40]=armttbarGL[14] + armttbarGL[47] - 145./16. + 
   armttbarGL[40];
   armttbarGL[42]=3*armttbarGL[13];
   armttbarGL[40]=armttbarGL[42]*armttbarGL[40];
   armttbarGL[45]=armttbarGL[10]*armttbarGL[13];
   armttbarGL[46]=9*armttbarGL[13];
   armttbarGL[46]=armttbarGL[15]*armttbarGL[46];
   armttbarGL[40]=armttbarGL[46] + 27*armttbarGL[45] + 3*armttbarGL[7]
    + 5./4.*armttbarGL[9] + armttbarGL[40];
   armttbarGL[46]= - 69./2. + armttbarGL[16];
   armttbarGL[46]=1./4.*armttbarGL[46] + armttbarGL[14];
   armttbarGL[47]=203./8.*armttbarGL[13] + armttbarGL[12];
   armttbarGL[47]=armttbarGL[3]*armttbarGL[47];
   armttbarGL[46]=armttbarGL[47] + 3./4.*armttbarGL[46];
   armttbarGL[46]=armttbarGL[12]*armttbarGL[46];
   armttbarGL[39]=85./16.*armttbarGL[19] + 17./16.*armttbarGL[25] + 1./
   2.*armttbarGL[40] + armttbarGL[46] + armttbarGL[39];
   armttbarGL[40]=pow(MMH,2);
   armttbarGL[46]=MMt*MMH;
   armttbarGL[46]= - 3./8.*armttbarGL[40] + armttbarGL[46];
   armttbarGL[46]=armttbarGL[22]*armttbarGL[46];
   armttbarGL[39]=3./4.*armttbarGL[46] + 1./4.*armttbarGL[39];
   armttbarGL[38]=armttbarGL[38]*armttbarGL[41];
   armttbarGL[39]=armttbarGL[38]*armttbarGL[39];
   armttbarGL[37]=armttbarGL[37] + armttbarGL[39];
   armttbarGL[37]=MMt*armttbarGL[37];
   armttbarGL[39]= - 13 + armttbarGL[24];
   armttbarGL[41]= - 15 + armttbarGL[11];
   armttbarGL[41]=armttbarGL[11]*armttbarGL[41];
   armttbarGL[46]= - armttbarGL[5]*armttbarGL[8];
   armttbarGL[47]= - armttbarGL[36] - 1 - armttbarGL[17];
   armttbarGL[47]=MMH*armttbarGL[5]*armttbarGL[47];
   armttbarGL[39]=armttbarGL[47] + armttbarGL[35] + armttbarGL[46] + 3./
   2.*armttbarGL[39] + armttbarGL[41];
   armttbarGL[41]= - 3./2.*armttbarGL[15] - 9./2.*armttbarGL[10];
   armttbarGL[41]=armttbarGL[43]*armttbarGL[41];
   armttbarGL[43]=1./2.*armttbarGL[5];
   armttbarGL[46]=armttbarGL[43]*armttbarGL[13];
   armttbarGL[47]=armttbarGL[20]*armttbarGL[43];
   armttbarGL[39]=armttbarGL[46] + armttbarGL[47] + 5./4.*
   armttbarGL[17] - 23./4.*armttbarGL[21] + 3*armttbarGL[26] + 
   armttbarGL[41] + 1./2.*armttbarGL[39];
   armttbarGL[39]=MMH*armttbarGL[39];
   armttbarGL[41]= - armttbarGL[5]*armttbarGL[42];
   armttbarGL[41]=armttbarGL[41] + 31 - 23*armttbarGL[11];
   armttbarGL[41]=armttbarGL[13]*armttbarGL[41];
   armttbarGL[42]=armttbarGL[15]*armttbarGL[42];
   armttbarGL[41]=armttbarGL[42] + 9*armttbarGL[45] + armttbarGL[41] + 
   31*armttbarGL[20] - 5./2.*armttbarGL[8] - 13*armttbarGL[7] - 3./2.*
   armttbarGL[6];
   armttbarGL[39]=1./2.*armttbarGL[41] + armttbarGL[39];
   armttbarGL[39]=MMH*armttbarGL[39];
   armttbarGL[41]=MMH*armttbarGL[43];
   armttbarGL[41]=armttbarGL[41] + 3./4.*armttbarGL[15] + 9./4.*
   armttbarGL[10] - armttbarGL[46] + 9 - 1./2.*armttbarGL[11];
   armttbarGL[41]=MMH*armttbarGL[41];
   armttbarGL[41]=1./8.*armttbarGL[12] - 51./2.*armttbarGL[13] + 
   armttbarGL[41];
   armttbarGL[41]=armttbarGL[12]*armttbarGL[41];
   armttbarGL[39]=armttbarGL[39] + armttbarGL[41];
   armttbarGL[41]=3./32.*armttbarGL[34] - 21./64.*armttbarGL[18];
   armttbarGL[40]=armttbarGL[40]*armttbarGL[41];
   armttbarGL[41]=armttbarGL[19]*MMH;
   armttbarGL[39]=17./64.*armttbarGL[41] + armttbarGL[40] + 
   armttbarGL[44] + 1./16.*armttbarGL[39];
   armttbarGL[38]=armttbarGL[39]*armttbarGL[38];
   armttbarGL[37]=armttbarGL[38] + armttbarGL[37];

      mttbarGLret = armttbarGL[37]*armttbarGL[1];
      return mttbarGLret.real();
}
