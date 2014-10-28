#include <tt.hpp>
std::complex<long double>
tt::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbar[61], mttbarret;

    armttbar[1]=double(nH);
    armttbar[2]=pow(CW,-1);
    armttbar[3]=pow(MMH,-1);
    armttbar[4]=pow(MMZ,-1);
    armttbar[5]=pow(SW,-1);
    armttbar[6]=Tsil::A(MMt,mu2);
    armttbar[7]=double(boson);
    armttbar[8]=pow(MMt,-1);
    armttbar[9]=Tsil::I2(MMH,MMt,MMt,mu2);
    armttbar[10]=Tsil::I2(MMZ,MMt,MMt,mu2);
    armttbar[11]=Tsil::I2(0,MMW,MMt,mu2);
    armttbar[12]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbar[13]=Tsil::A(MMH,mu2);
    armttbar[14]=Tsil::B(MMZ,MMt,MMt,mu2);
    armttbar[15]=Tsil::A(MMZ,mu2);
    armttbar[16]=Tsil::Beps(MMH,MMt,MMt,mu2);
    armttbar[17]=Tsil::Beps(MMZ,MMt,MMt,mu2);
    armttbar[18]=Tsil::A(MMW,mu2);
    armttbar[19]=std::real(Tsil::B(0,MMW,MMt,mu2));
    armttbar[20]=Tsil::Aeps(MMH,mu2);
    armttbar[21]=Tsil::Aeps(MMZ,mu2);
    armttbar[22]=Tsil::Aeps(MMW,mu2);
    armttbar[23]=Tsil::Aeps(MMt,mu2);
    armttbar[24]=std::real(Tsil::Beps(0,MMW,MMt,mu2));
    armttbar[25]=protWt000->M(0);
    armttbar[26]=prot0ttHt->M(0);
    armttbar[27]=prot0ttZt->M(0);
    armttbar[28]=prot0tt0t->M(0);
    armttbar[29]=prottH0H->Vxzuv(0);
    armttbar[30]=prottZ0Z->Vxzuv(0);
    armttbar[31]=protWt000->Uzxyv(0);
    armttbar[32]=prot0ttHt->Tuxv(0);
    armttbar[33]=prot0ttZt->Tuxv(0);
    armttbar[34]=protWt000->Txuv(0);
    armttbar[35]=protWt000->Tyzv(0);
    armttbar[36]=1/(4*MMt - MMZ);
   armttbar[37]=pow(Pi,2);
   armttbar[38]= - 5./2.*armttbar[37];
   armttbar[39]=5./2.*armttbar[34] + armttbar[38] + 13*armttbar[32] - 
   29./4. - 8*armttbar[16];
   armttbar[40]=1 + armttbar[19];
   armttbar[40]=armttbar[19]*armttbar[40];
   armttbar[41]=armttbar[3]*armttbar[20];
   armttbar[42]= - 2*armttbar[14];
   armttbar[43]= - 4./3.*armttbar[29] + armttbar[26];
   armttbar[43]=MMH*armttbar[43];
   armttbar[44]= - 8*armttbar[26] - armttbar[25];
   armttbar[44]=MMt*armttbar[44];
   armttbar[39]=2./3.*armttbar[44] + 4*armttbar[43] + armttbar[42] + 13.
   /3.*armttbar[41] + 5./6.*armttbar[40] + 52./3.*armttbar[12] - 
   armttbar[33] + 1./3.*armttbar[39] + 3*armttbar[35];
   armttbar[39]=MMt*armttbar[39];
   armttbar[40]= - armttbar[3]*armttbar[13];
   armttbar[41]=4./3.*armttbar[19];
   armttbar[40]=13./3.*armttbar[40] + armttbar[41] + 23./2. + 14./3.*
   armttbar[12];
   armttbar[43]=armttbar[15]*armttbar[36];
   armttbar[44]=100./9.*armttbar[43];
   armttbar[45]= - armttbar[6]*armttbar[36];
   armttbar[46]=100./3.*armttbar[45] + armttbar[40] + armttbar[44];
   armttbar[46]=armttbar[6]*armttbar[46];
   armttbar[47]= - 13 + 4*armttbar[16];
   armttbar[48]=2*armttbar[29] - armttbar[26];
   armttbar[48]=MMH*armttbar[48];
   armttbar[47]=2./3.*armttbar[48] - 17./3.*armttbar[12] - 4./3.*
   armttbar[35] + 1./3.*armttbar[47] - 2*armttbar[32];
   armttbar[47]=MMH*armttbar[47];
   armttbar[44]= - 1./2. + armttbar[44];
   armttbar[44]=armttbar[15]*armttbar[44];
   armttbar[48]= - 13./6.*armttbar[11];
   armttbar[49]=7./2.*armttbar[13] + 17./6.*armttbar[22] + armttbar[48]
    + 3*armttbar[20] - 11./3.*armttbar[9];
   armttbar[50]= - 2./3.*armttbar[12]*armttbar[13];
   armttbar[51]=59./6.*armttbar[23];
   armttbar[52]= - 1./3. + 1./2.*armttbar[19];
   armttbar[52]=5*armttbar[18]*armttbar[52];
   armttbar[53]= - 4./3.*armttbar[3]*pow(armttbar[13],2);
   armttbar[44]=armttbar[46] + armttbar[39] + armttbar[44] + 
   armttbar[47] + armttbar[53] + armttbar[52] + armttbar[51] + 
   armttbar[50] + armttbar[49] + 65./27.*armttbar[21];
   armttbar[46]=pow(armttbar[2],2);
   armttbar[44]=armttbar[46]*armttbar[44];
   armttbar[43]=4*armttbar[43];
   armttbar[40]=12*armttbar[45] + armttbar[40] + armttbar[43];
   armttbar[40]=armttbar[6]*armttbar[40];
   armttbar[43]= - 1./2. + armttbar[43];
   armttbar[43]=armttbar[15]*armttbar[43];
   armttbar[39]=armttbar[40] + armttbar[39] + armttbar[43] + 
   armttbar[47] + armttbar[53] + armttbar[52] + armttbar[51] + 
   armttbar[50] + armttbar[49] + 1./3.*armttbar[21];
   armttbar[40]=pow(armttbar[5],2);
   armttbar[39]=armttbar[40]*armttbar[39];
   armttbar[43]= - 1 + armttbar[19];
   armttbar[43]=armttbar[18]*armttbar[43];
   armttbar[47]=5./6.*armttbar[43] + 5./6.*armttbar[22] - armttbar[21];
   armttbar[47]=MMt*armttbar[47];
   armttbar[49]=armttbar[6]*armttbar[15];
   armttbar[50]=pow(armttbar[18],2);
   armttbar[47]=armttbar[49] + 5./6.*armttbar[50] + armttbar[47];
   armttbar[43]=armttbar[22] + armttbar[43];
   armttbar[43]=MMt*armttbar[43];
   armttbar[43]=armttbar[50] + armttbar[43];
   armttbar[43]=armttbar[46]*armttbar[43];
   armttbar[43]=armttbar[47] + 5./6.*armttbar[43];
   armttbar[43]=armttbar[46]*armttbar[43];
   armttbar[47]=armttbar[40]*armttbar[47];
   armttbar[43]=armttbar[43] + armttbar[47];
   armttbar[43]=armttbar[4]*armttbar[43];
   armttbar[47]= - 7./3.*armttbar[21];
   armttbar[51]=pow(armttbar[15],2);
   armttbar[52]= - armttbar[51]*armttbar[36];
   armttbar[52]=armttbar[47] + 8*armttbar[52];
   armttbar[53]= - armttbar[15]*armttbar[36];
   armttbar[54]=armttbar[6]*armttbar[36];
   armttbar[53]=1./3.*armttbar[53] + armttbar[54];
   armttbar[53]=armttbar[6]*armttbar[53];
   armttbar[52]=1./3.*armttbar[52] + 8*armttbar[53];
   armttbar[39]=armttbar[43] + armttbar[39] + 32./3.*armttbar[52] + 
   armttbar[44];
   armttbar[39]=armttbar[4]*armttbar[39];
   armttbar[43]= - armttbar[3] + armttbar[27] - 2./3.*armttbar[30] + 
   armttbar[25];
   armttbar[44]=2*armttbar[14];
   armttbar[52]=armttbar[44] + 13./2. - 2*armttbar[17];
   armttbar[52]=armttbar[36]*armttbar[52];
   armttbar[43]=2*armttbar[43] + armttbar[52];
   armttbar[43]=MMZ*armttbar[43];
   armttbar[53]= - 119./8. + 2*armttbar[34];
   armttbar[55]=5*armttbar[19];
   armttbar[56]= - 1 + armttbar[55];
   armttbar[56]=armttbar[19]*armttbar[56];
   armttbar[57]= - armttbar[3]*armttbar[18];
   armttbar[58]=2*armttbar[23] + 3*armttbar[10] - 7*armttbar[21];
   armttbar[58]=armttbar[36]*armttbar[58];
   armttbar[44]= - 5 + armttbar[44];
   armttbar[44]=armttbar[36]*armttbar[44];
   armttbar[59]= - armttbar[3] + armttbar[44];
   armttbar[59]=armttbar[15]*armttbar[59];
   armttbar[60]= - MMt*armttbar[27];
   armttbar[43]=armttbar[43] + 6*armttbar[45] + 4./3.*armttbar[60] + 
   armttbar[59] + armttbar[58] + 10./3.*armttbar[14] + 2*armttbar[57]
    + 1./3.*armttbar[56] + 1./3.*armttbar[33] - 2./3.*armttbar[17] + 1./
   3.*armttbar[53] + 2*armttbar[35];
   armttbar[43]=armttbar[40]*armttbar[43];
   armttbar[53]=armttbar[3]*armttbar[18];
   armttbar[56]=3*armttbar[15]*armttbar[3];
   armttbar[41]=armttbar[56] + 7./6.*armttbar[14] + 6*armttbar[53] + 3./
   2. + armttbar[41];
   armttbar[41]=armttbar[6]*armttbar[41];
   armttbar[53]=armttbar[6]*armttbar[3];
   armttbar[57]=armttbar[30] - armttbar[25];
   armttbar[57]=2*armttbar[57] - armttbar[27];
   armttbar[57]=MMZ*armttbar[57];
   armttbar[57]=2./3.*armttbar[57] + 6*armttbar[53] - 7./6.*
   armttbar[14] - 47./6.*armttbar[19] + 3./2.*armttbar[33] + 4./3.*
   armttbar[17] - 19./6.*armttbar[34] + 5./2.*armttbar[37] - 2*
   armttbar[24] - 1 + 2*armttbar[31];
   armttbar[57]=MMZ*armttbar[57];
   armttbar[55]=11./2. + armttbar[55];
   armttbar[55]=armttbar[18]*armttbar[55];
   armttbar[58]=5./2. - 2./3.*armttbar[14];
   armttbar[58]=armttbar[15]*armttbar[58];
   armttbar[41]=armttbar[57] + armttbar[41] + armttbar[58] + 1./3.*
   armttbar[55] + 23./6.*armttbar[23] + 4*armttbar[21] - 7./6.*
   armttbar[10] + armttbar[48] - armttbar[22];
   armttbar[41]=armttbar[40]*armttbar[41];
   armttbar[48]=pow(CW,2);
   armttbar[55]=1 + armttbar[48];
   armttbar[57]=armttbar[37]*armttbar[55];
   armttbar[58]= - 1 - armttbar[48];
   armttbar[59]=armttbar[34]*armttbar[58];
   armttbar[58]=armttbar[19]*armttbar[58];
   armttbar[48]=32./9.*armttbar[14] + armttbar[58] + 32./9.*
   armttbar[33] + armttbar[59] + 5*armttbar[57] + 23./9. - armttbar[48]
   ;
   armttbar[48]=MMZ*armttbar[48];
   armttbar[57]=56./9.*armttbar[14] + 22./9. + armttbar[19];
   armttbar[57]=armttbar[6]*armttbar[57];
   armttbar[48]=armttbar[48] + 2*armttbar[57] - 32./9.*armttbar[15] + 
   armttbar[18] + 37./9.*armttbar[23] - armttbar[11] + 32./9.*
   armttbar[10];
   armttbar[48]=MMZ*armttbar[48];
   armttbar[57]= - 7./2.*armttbar[14];
   armttbar[58]=2 + armttbar[57];
   armttbar[58]=armttbar[6]*armttbar[58];
   armttbar[59]= - armttbar[14] - 1 - armttbar[33];
   armttbar[59]=MMZ*armttbar[59];
   armttbar[58]=armttbar[59] + armttbar[58] + armttbar[15] - 
   armttbar[10] - 2*armttbar[23];
   armttbar[58]=MMZ*armttbar[58];
   armttbar[59]= - 85./2.*armttbar[15] + 41*armttbar[6];
   armttbar[59]=armttbar[6]*armttbar[59];
   armttbar[58]=armttbar[59] + 17./3.*armttbar[58];
   armttbar[58]=armttbar[46]*armttbar[58];
   armttbar[48]=1./3.*armttbar[58] + 80./3.*armttbar[49] + armttbar[48]
   ;
   armttbar[49]=pow(armttbar[6],2);
   armttbar[58]=armttbar[49]*MMH;
   armttbar[59]=armttbar[46]*armttbar[58];
   armttbar[58]=armttbar[40]*armttbar[58];
   armttbar[58]=armttbar[59] + armttbar[58];
   armttbar[58]=armttbar[4]*armttbar[58];
   armttbar[59]=MMZ*armttbar[49];
   armttbar[60]=armttbar[46]*armttbar[59];
   armttbar[49]= - MMZ*armttbar[49];
   armttbar[49]=32*armttbar[49] + 17*armttbar[60];
   armttbar[59]=armttbar[40]*armttbar[59];
   armttbar[49]=1./9.*armttbar[49] + armttbar[59];
   armttbar[49]=armttbar[8]*armttbar[49];
   armttbar[59]= - 5 - armttbar[19];
   armttbar[57]=2*armttbar[59] + armttbar[57];
   armttbar[57]=armttbar[6]*armttbar[57];
   armttbar[37]= - armttbar[14] + armttbar[19] - armttbar[33] - 5*
   armttbar[37] + armttbar[34];
   armttbar[37]=MMZ*armttbar[37];
   armttbar[37]=armttbar[37] + armttbar[57] + armttbar[15] - 
   armttbar[18] + armttbar[23] + armttbar[11] - armttbar[10];
   armttbar[37]=MMZ*armttbar[37];
   armttbar[57]= - 1./3.*armttbar[18] - 1./2.*armttbar[15];
   armttbar[57]=5*armttbar[57] + armttbar[6];
   armttbar[57]=armttbar[6]*armttbar[57];
   armttbar[37]=armttbar[57] + 1./3.*armttbar[37];
   armttbar[37]=armttbar[40]*armttbar[37];
   armttbar[37]=armttbar[49] + armttbar[58] + 1./3.*armttbar[48] + 
   armttbar[37];
   armttbar[37]=armttbar[8]*armttbar[37];
   armttbar[48]=13./9.*armttbar[15] - 2*armttbar[6];
   armttbar[48]=armttbar[6]*armttbar[48];
   armttbar[48]=4./9.*armttbar[51] + armttbar[48];
   armttbar[49]=armttbar[9] - armttbar[13];
   armttbar[57]=armttbar[12] + 1 + armttbar[32];
   armttbar[57]=MMH*armttbar[57];
   armttbar[49]=1./2.*armttbar[57] + 1./2.*armttbar[49] - 5*
   armttbar[23];
   armttbar[49]=MMH*armttbar[49];
   armttbar[49]=5*armttbar[50] + armttbar[49];
   armttbar[50]=armttbar[49] - 68./9.*armttbar[51];
   armttbar[57]= - 1 - 7./2.*armttbar[12];
   armttbar[57]=MMH*armttbar[57];
   armttbar[57]= - 3./2.*armttbar[18] + 1./3.*armttbar[57];
   armttbar[58]=19./3.*armttbar[6] + armttbar[57] - 140./27.*
   armttbar[15];
   armttbar[58]=armttbar[6]*armttbar[58];
   armttbar[50]=1./3.*armttbar[50] + armttbar[58];
   armttbar[50]=armttbar[46]*armttbar[50];
   armttbar[49]=armttbar[49] - 4*armttbar[51];
   armttbar[51]=armttbar[6] + armttbar[57] - 4./3.*armttbar[15];
   armttbar[51]=armttbar[6]*armttbar[51];
   armttbar[49]=1./3.*armttbar[49] + armttbar[51];
   armttbar[49]=armttbar[40]*armttbar[49];
   armttbar[48]=armttbar[49] + 32./3.*armttbar[48] + armttbar[50];
   armttbar[48]=armttbar[4]*armttbar[48];
   armttbar[49]=1./3.*armttbar[10] - 2*armttbar[21];
   armttbar[50]= - 1 + 1./3.*armttbar[14];
   armttbar[50]=armttbar[15]*armttbar[50];
   armttbar[51]=10 + 7*armttbar[14];
   armttbar[51]=armttbar[6]*armttbar[51];
   armttbar[49]=1./3.*armttbar[51] + 2*armttbar[50] + 2*armttbar[49] - 
   armttbar[23];
   armttbar[50]= - armttbar[6]*armttbar[3];
   armttbar[51]=armttbar[25]*armttbar[55];
   armttbar[51]=16./9.*armttbar[27] - 32./9.*armttbar[30] + 
   armttbar[51];
   armttbar[51]=MMZ*armttbar[51];
   armttbar[38]=4./3.*armttbar[51] + 4*armttbar[50] + 64./27.*
   armttbar[14] + 47./6.*armttbar[19] - 64./9.*armttbar[33] - 128./27.*
   armttbar[17] + 92./27.*armttbar[35] + 19./6.*armttbar[34] + 
   armttbar[38] + 2*armttbar[24] + 211./54. - 2*armttbar[31];
   armttbar[38]=MMZ*armttbar[38];
   armttbar[50]=68*armttbar[17] - 163./2. - 68*armttbar[35];
   armttbar[50]= - 95./6.*armttbar[14] + 1./3.*armttbar[50] + 59./2.*
   armttbar[33];
   armttbar[51]=2*armttbar[30] - armttbar[27];
   armttbar[51]=MMZ*armttbar[51];
   armttbar[50]=34./27.*armttbar[51] + 1./9.*armttbar[50] + 2*
   armttbar[53];
   armttbar[50]=MMZ*armttbar[50];
   armttbar[51]=17 - 7*armttbar[14];
   armttbar[51]=7./54.*armttbar[51] + armttbar[56];
   armttbar[51]=armttbar[6]*armttbar[51];
   armttbar[53]=77./2. - 34./3.*armttbar[14];
   armttbar[53]=armttbar[15]*armttbar[53];
   armttbar[53]=armttbar[53] - 23*armttbar[23] - 95./6.*armttbar[10] + 
   68*armttbar[21];
   armttbar[50]=armttbar[50] + 1./9.*armttbar[53] + armttbar[51];
   armttbar[50]=armttbar[46]*armttbar[50];
   armttbar[37]=armttbar[37] + armttbar[48] + armttbar[41] + 
   armttbar[50] + 32./9.*armttbar[49] + armttbar[38];
   armttbar[37]=armttbar[8]*armttbar[37];
   armttbar[38]=11*armttbar[14] - 7*armttbar[33] - 4*armttbar[17] - 22
    - 7*armttbar[35];
   armttbar[41]= - 2./3.*armttbar[23] - armttbar[10] + 7./3.*
   armttbar[21];
   armttbar[41]=armttbar[36]*armttbar[41];
   armttbar[42]=5 + armttbar[42];
   armttbar[42]=armttbar[15]*armttbar[36]*armttbar[42];
   armttbar[48]=armttbar[28] - armttbar[27];
   armttbar[48]=MMt*armttbar[48];
   armttbar[38]=4*armttbar[54] + 8./9.*armttbar[48] + 2./3.*
   armttbar[42] + 1./9.*armttbar[38] + 2*armttbar[41];
   armttbar[41]= - 4*armttbar[14] - 13 + 4*armttbar[17];
   armttbar[41]=armttbar[36]*armttbar[41];
   armttbar[41]=16./9.*armttbar[41] + 2./3.*armttbar[3] - 128./27.*
   armttbar[30] - armttbar[25];
   armttbar[41]=MMZ*armttbar[41];
   armttbar[38]=16./3.*armttbar[38] + armttbar[41];
   armttbar[41]=2./3.*armttbar[14] + 65./3.*armttbar[33] + 14./3.*
   armttbar[17] - 47./8. - 28./3.*armttbar[35];
   armttbar[42]=2./3.*armttbar[23] + armttbar[10] + armttbar[47];
   armttbar[42]=armttbar[36]*armttbar[42];
   armttbar[41]=1./3.*armttbar[41] + 25*armttbar[42];
   armttbar[42]= - 1./3.*armttbar[3] + 14./27.*armttbar[30] + 
   armttbar[27];
   armttbar[42]=2*armttbar[42] + 25./9.*armttbar[52];
   armttbar[42]=MMZ*armttbar[42];
   armttbar[44]= - armttbar[3] + 25./9.*armttbar[44];
   armttbar[44]=armttbar[15]*armttbar[44];
   armttbar[47]=MMt*armttbar[27];
   armttbar[41]=armttbar[42] + 50./3.*armttbar[45] + 28./27.*
   armttbar[47] + 1./3.*armttbar[41] + armttbar[44];
   armttbar[41]=armttbar[46]*armttbar[41];
   armttbar[37]=armttbar[37] + armttbar[39] + armttbar[43] + 2*
   armttbar[38] + armttbar[41];
   armttbar[37]=armttbar[7]*armttbar[37];
   armttbar[38]=armttbar[3]*armttbar[1];
   armttbar[39]=MMt*armttbar[38];
   armttbar[41]= - armttbar[6]*armttbar[3]*armttbar[1];
   armttbar[39]=armttbar[39] + 9*armttbar[41];
   armttbar[39]=armttbar[6]*armttbar[39];
   armttbar[38]=armttbar[38]*pow(MMt,2);
   armttbar[38]=8*armttbar[38] + armttbar[39];
   armttbar[39]=armttbar[46]*armttbar[38];
   armttbar[38]=armttbar[40]*armttbar[38];
   armttbar[38]=armttbar[39] + armttbar[38];
   armttbar[38]=armttbar[4]*armttbar[38];

      mttbarret = armttbar[37] + 4*armttbar[38];
      return mttbarret;
}
