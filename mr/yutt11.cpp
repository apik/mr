#include <tt.hpp>
std::complex<long double>
tt::my11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryutt[61], yuttret;

    aryutt[1]=double(nH);
    aryutt[2]=pow(CW,-1);
    aryutt[3]=pow(MMH,-1);
    aryutt[4]=pow(MMZ,-1);
    aryutt[5]=pow(SW,-1);
    aryutt[6]=Tsil::A(MMt,mu2);
    aryutt[7]=double(boson);
    aryutt[8]=pow(MMt,-1);
    aryutt[9]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryutt[10]=Tsil::I2(MMZ,MMt,MMt,mu2);
    aryutt[11]=Tsil::I2(0,MMW,MMt,mu2);
    aryutt[12]=Tsil::B(MMH,MMt,MMt,mu2);
    aryutt[13]=Tsil::A(MMH,mu2);
    aryutt[14]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryutt[15]=Tsil::A(MMZ,mu2);
    aryutt[16]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryutt[17]=Tsil::Beps(MMZ,MMt,MMt,mu2);
    aryutt[18]=Tsil::A(MMW,mu2);
    aryutt[19]=std::real(Tsil::B(0,MMW,MMt,mu2));
    aryutt[20]=Tsil::Aeps(MMH,mu2);
    aryutt[21]=Tsil::Aeps(MMZ,mu2);
    aryutt[22]=Tsil::Aeps(MMW,mu2);
    aryutt[23]=Tsil::Aeps(MMt,mu2);
    aryutt[24]=std::real(Tsil::Beps(0,MMW,MMt,mu2));
    aryutt[25]=protWt000->M(0);
    aryutt[26]=prot0ttHt->M(0);
    aryutt[27]=prot0ttZt->M(0);
    aryutt[28]=prot0tt0t->M(0);
    aryutt[29]=prottH0H->Vxzuv(0);
    aryutt[30]=prottZ0Z->Vxzuv(0);
    aryutt[31]=protWt000->Uzxyv(0);
    aryutt[32]=prot0ttHt->Tuxv(0);
    aryutt[33]=prot0ttZt->Tuxv(0);
    aryutt[34]=protWt000->Txuv(0);
    aryutt[35]=protWt000->Tyzv(0);
    aryutt[36]=1/(4*MMt - MMZ);
   aryutt[37]=pow(Pi,2);
   aryutt[38]= - 5./2.*aryutt[37];
   aryutt[39]=5./2.*aryutt[34] + aryutt[38] + 13*aryutt[32] - 29./4. - 
   8*aryutt[16];
   aryutt[40]=1 + aryutt[19];
   aryutt[40]=aryutt[19]*aryutt[40];
   aryutt[41]=aryutt[3]*aryutt[20];
   aryutt[42]= - 2*aryutt[14];
   aryutt[43]= - 4./3.*aryutt[29] + aryutt[26];
   aryutt[43]=MMH*aryutt[43];
   aryutt[44]= - 8*aryutt[26] - aryutt[25];
   aryutt[44]=MMt*aryutt[44];
   aryutt[39]=2./3.*aryutt[44] + 4*aryutt[43] + aryutt[42] + 13./3.*
   aryutt[41] + 5./6.*aryutt[40] + 52./3.*aryutt[12] - aryutt[33] + 1./
   3.*aryutt[39] + 3*aryutt[35];
   aryutt[39]=MMt*aryutt[39];
   aryutt[40]= - aryutt[3]*aryutt[13];
   aryutt[41]=4./3.*aryutt[19];
   aryutt[40]=13./3.*aryutt[40] + aryutt[41] + 23./2. + 14./3.*
   aryutt[12];
   aryutt[43]=aryutt[15]*aryutt[36];
   aryutt[44]=100./9.*aryutt[43];
   aryutt[45]= - aryutt[6]*aryutt[36];
   aryutt[46]=100./3.*aryutt[45] + aryutt[40] + aryutt[44];
   aryutt[46]=aryutt[6]*aryutt[46];
   aryutt[47]= - 13 + 4*aryutt[16];
   aryutt[48]=2*aryutt[29] - aryutt[26];
   aryutt[48]=MMH*aryutt[48];
   aryutt[47]=2./3.*aryutt[48] - 17./3.*aryutt[12] - 4./3.*aryutt[35]
    + 1./3.*aryutt[47] - 2*aryutt[32];
   aryutt[47]=MMH*aryutt[47];
   aryutt[44]= - 1./2. + aryutt[44];
   aryutt[44]=aryutt[15]*aryutt[44];
   aryutt[48]= - 13./6.*aryutt[11];
   aryutt[49]=7./2.*aryutt[13] + 17./6.*aryutt[22] + aryutt[48] + 3*
   aryutt[20] - 11./3.*aryutt[9];
   aryutt[50]= - 2./3.*aryutt[12]*aryutt[13];
   aryutt[51]=59./6.*aryutt[23];
   aryutt[52]= - 1./3. + 1./2.*aryutt[19];
   aryutt[52]=5*aryutt[18]*aryutt[52];
   aryutt[53]= - 4./3.*aryutt[3]*pow(aryutt[13],2);
   aryutt[44]=aryutt[46] + aryutt[39] + aryutt[44] + aryutt[47] + 
   aryutt[53] + aryutt[52] + aryutt[51] + aryutt[50] + aryutt[49] + 65./
   27.*aryutt[21];
   aryutt[46]=pow(aryutt[2],2);
   aryutt[44]=aryutt[46]*aryutt[44];
   aryutt[43]=4*aryutt[43];
   aryutt[40]=12*aryutt[45] + aryutt[40] + aryutt[43];
   aryutt[40]=aryutt[6]*aryutt[40];
   aryutt[43]= - 1./2. + aryutt[43];
   aryutt[43]=aryutt[15]*aryutt[43];
   aryutt[39]=aryutt[40] + aryutt[39] + aryutt[43] + aryutt[47] + 
   aryutt[53] + aryutt[52] + aryutt[51] + aryutt[50] + aryutt[49] + 1./
   3.*aryutt[21];
   aryutt[40]=pow(aryutt[5],2);
   aryutt[39]=aryutt[40]*aryutt[39];
   aryutt[43]= - 1 + aryutt[19];
   aryutt[43]=aryutt[18]*aryutt[43];
   aryutt[47]=5./6.*aryutt[43] + 5./6.*aryutt[22] - aryutt[21];
   aryutt[47]=MMt*aryutt[47];
   aryutt[49]=aryutt[6]*aryutt[15];
   aryutt[50]=pow(aryutt[18],2);
   aryutt[47]=aryutt[49] + 5./6.*aryutt[50] + aryutt[47];
   aryutt[43]=aryutt[22] + aryutt[43];
   aryutt[43]=MMt*aryutt[43];
   aryutt[43]=aryutt[50] + aryutt[43];
   aryutt[43]=aryutt[46]*aryutt[43];
   aryutt[43]=aryutt[47] + 5./6.*aryutt[43];
   aryutt[43]=aryutt[46]*aryutt[43];
   aryutt[47]=aryutt[40]*aryutt[47];
   aryutt[43]=aryutt[43] + aryutt[47];
   aryutt[43]=aryutt[4]*aryutt[43];
   aryutt[47]= - 7./3.*aryutt[21];
   aryutt[51]=pow(aryutt[15],2);
   aryutt[52]= - aryutt[51]*aryutt[36];
   aryutt[52]=aryutt[47] + 8*aryutt[52];
   aryutt[53]= - aryutt[15]*aryutt[36];
   aryutt[54]=aryutt[6]*aryutt[36];
   aryutt[53]=1./3.*aryutt[53] + aryutt[54];
   aryutt[53]=aryutt[6]*aryutt[53];
   aryutt[52]=1./3.*aryutt[52] + 8*aryutt[53];
   aryutt[39]=aryutt[43] + aryutt[39] + 32./3.*aryutt[52] + aryutt[44];
   aryutt[39]=aryutt[4]*aryutt[39];
   aryutt[43]= - aryutt[3] + aryutt[27] - 2./3.*aryutt[30] + aryutt[25]
   ;
   aryutt[44]=2*aryutt[14];
   aryutt[52]=aryutt[44] + 13./2. - 2*aryutt[17];
   aryutt[52]=aryutt[36]*aryutt[52];
   aryutt[43]=2*aryutt[43] + aryutt[52];
   aryutt[43]=MMZ*aryutt[43];
   aryutt[53]= - 119./8. + 2*aryutt[34];
   aryutt[55]=5*aryutt[19];
   aryutt[56]= - 1 + aryutt[55];
   aryutt[56]=aryutt[19]*aryutt[56];
   aryutt[57]= - aryutt[3]*aryutt[18];
   aryutt[58]=2*aryutt[23] + 3*aryutt[10] - 7*aryutt[21];
   aryutt[58]=aryutt[36]*aryutt[58];
   aryutt[44]= - 5 + aryutt[44];
   aryutt[44]=aryutt[36]*aryutt[44];
   aryutt[59]= - aryutt[3] + aryutt[44];
   aryutt[59]=aryutt[15]*aryutt[59];
   aryutt[60]= - MMt*aryutt[27];
   aryutt[43]=aryutt[43] + 6*aryutt[45] + 4./3.*aryutt[60] + aryutt[59]
    + aryutt[58] + 10./3.*aryutt[14] + 2*aryutt[57] + 1./3.*aryutt[56]
    + 1./3.*aryutt[33] - 2./3.*aryutt[17] + 1./3.*aryutt[53] + 2*
   aryutt[35];
   aryutt[43]=aryutt[40]*aryutt[43];
   aryutt[53]=aryutt[3]*aryutt[18];
   aryutt[56]=3*aryutt[15]*aryutt[3];
   aryutt[41]=aryutt[56] + 7./6.*aryutt[14] + 6*aryutt[53] + 3./2. + 
   aryutt[41];
   aryutt[41]=aryutt[6]*aryutt[41];
   aryutt[53]=aryutt[6]*aryutt[3];
   aryutt[57]=aryutt[30] - aryutt[25];
   aryutt[57]=2*aryutt[57] - aryutt[27];
   aryutt[57]=MMZ*aryutt[57];
   aryutt[57]=2./3.*aryutt[57] + 6*aryutt[53] - 7./6.*aryutt[14] - 47./
   6.*aryutt[19] + 3./2.*aryutt[33] + 4./3.*aryutt[17] - 19./6.*
   aryutt[34] + 5./2.*aryutt[37] - 2*aryutt[24] - 1 + 2*aryutt[31];
   aryutt[57]=MMZ*aryutt[57];
   aryutt[55]=11./2. + aryutt[55];
   aryutt[55]=aryutt[18]*aryutt[55];
   aryutt[58]=5./2. - 2./3.*aryutt[14];
   aryutt[58]=aryutt[15]*aryutt[58];
   aryutt[41]=aryutt[57] + aryutt[41] + aryutt[58] + 1./3.*aryutt[55]
    + 23./6.*aryutt[23] + 4*aryutt[21] - 7./6.*aryutt[10] + aryutt[48]
    - aryutt[22];
   aryutt[41]=aryutt[40]*aryutt[41];
   aryutt[48]=pow(CW,2);
   aryutt[55]=1 + aryutt[48];
   aryutt[57]=aryutt[37]*aryutt[55];
   aryutt[58]= - 1 - aryutt[48];
   aryutt[59]=aryutt[34]*aryutt[58];
   aryutt[58]=aryutt[19]*aryutt[58];
   aryutt[48]=32./9.*aryutt[14] + aryutt[58] + 32./9.*aryutt[33] + 
   aryutt[59] + 5*aryutt[57] + 23./9. - aryutt[48];
   aryutt[48]=MMZ*aryutt[48];
   aryutt[57]=56./9.*aryutt[14] + 22./9. + aryutt[19];
   aryutt[57]=aryutt[6]*aryutt[57];
   aryutt[48]=aryutt[48] + 2*aryutt[57] - 32./9.*aryutt[15] + 
   aryutt[18] + 37./9.*aryutt[23] - aryutt[11] + 32./9.*aryutt[10];
   aryutt[48]=MMZ*aryutt[48];
   aryutt[57]= - 7./2.*aryutt[14];
   aryutt[58]=2 + aryutt[57];
   aryutt[58]=aryutt[6]*aryutt[58];
   aryutt[59]= - aryutt[14] - 1 - aryutt[33];
   aryutt[59]=MMZ*aryutt[59];
   aryutt[58]=aryutt[59] + aryutt[58] + aryutt[15] - aryutt[10] - 2*
   aryutt[23];
   aryutt[58]=MMZ*aryutt[58];
   aryutt[59]= - 85./2.*aryutt[15] + 41*aryutt[6];
   aryutt[59]=aryutt[6]*aryutt[59];
   aryutt[58]=aryutt[59] + 17./3.*aryutt[58];
   aryutt[58]=aryutt[46]*aryutt[58];
   aryutt[48]=1./3.*aryutt[58] + 80./3.*aryutt[49] + aryutt[48];
   aryutt[49]=pow(aryutt[6],2);
   aryutt[58]=aryutt[49]*MMH;
   aryutt[59]=aryutt[46]*aryutt[58];
   aryutt[58]=aryutt[40]*aryutt[58];
   aryutt[58]=aryutt[59] + aryutt[58];
   aryutt[58]=aryutt[4]*aryutt[58];
   aryutt[59]=MMZ*aryutt[49];
   aryutt[60]=aryutt[46]*aryutt[59];
   aryutt[49]= - MMZ*aryutt[49];
   aryutt[49]=32*aryutt[49] + 17*aryutt[60];
   aryutt[59]=aryutt[40]*aryutt[59];
   aryutt[49]=1./9.*aryutt[49] + aryutt[59];
   aryutt[49]=aryutt[8]*aryutt[49];
   aryutt[59]= - 5 - aryutt[19];
   aryutt[57]=2*aryutt[59] + aryutt[57];
   aryutt[57]=aryutt[6]*aryutt[57];
   aryutt[37]= - aryutt[14] + aryutt[19] - aryutt[33] - 5*aryutt[37] + 
   aryutt[34];
   aryutt[37]=MMZ*aryutt[37];
   aryutt[37]=aryutt[37] + aryutt[57] + aryutt[15] - aryutt[18] + 
   aryutt[23] + aryutt[11] - aryutt[10];
   aryutt[37]=MMZ*aryutt[37];
   aryutt[57]= - 1./3.*aryutt[18] - 1./2.*aryutt[15];
   aryutt[57]=5*aryutt[57] + aryutt[6];
   aryutt[57]=aryutt[6]*aryutt[57];
   aryutt[37]=aryutt[57] + 1./3.*aryutt[37];
   aryutt[37]=aryutt[40]*aryutt[37];
   aryutt[37]=aryutt[49] + aryutt[58] + 1./3.*aryutt[48] + aryutt[37];
   aryutt[37]=aryutt[8]*aryutt[37];
   aryutt[48]=13./9.*aryutt[15] - 2*aryutt[6];
   aryutt[48]=aryutt[6]*aryutt[48];
   aryutt[48]=4./9.*aryutt[51] + aryutt[48];
   aryutt[49]=aryutt[9] - aryutt[13];
   aryutt[57]=aryutt[12] + 1 + aryutt[32];
   aryutt[57]=MMH*aryutt[57];
   aryutt[49]=1./2.*aryutt[57] + 1./2.*aryutt[49] - 5*aryutt[23];
   aryutt[49]=MMH*aryutt[49];
   aryutt[49]=5*aryutt[50] + aryutt[49];
   aryutt[50]=aryutt[49] - 68./9.*aryutt[51];
   aryutt[57]= - 1 - 7./2.*aryutt[12];
   aryutt[57]=MMH*aryutt[57];
   aryutt[57]= - 3./2.*aryutt[18] + 1./3.*aryutt[57];
   aryutt[58]=19./3.*aryutt[6] + aryutt[57] - 140./27.*aryutt[15];
   aryutt[58]=aryutt[6]*aryutt[58];
   aryutt[50]=1./3.*aryutt[50] + aryutt[58];
   aryutt[50]=aryutt[46]*aryutt[50];
   aryutt[49]=aryutt[49] - 4*aryutt[51];
   aryutt[51]=aryutt[6] + aryutt[57] - 4./3.*aryutt[15];
   aryutt[51]=aryutt[6]*aryutt[51];
   aryutt[49]=1./3.*aryutt[49] + aryutt[51];
   aryutt[49]=aryutt[40]*aryutt[49];
   aryutt[48]=aryutt[49] + 32./3.*aryutt[48] + aryutt[50];
   aryutt[48]=aryutt[4]*aryutt[48];
   aryutt[49]=1./3.*aryutt[10] - 2*aryutt[21];
   aryutt[50]= - 1 + 1./3.*aryutt[14];
   aryutt[50]=aryutt[15]*aryutt[50];
   aryutt[51]=10 + 7*aryutt[14];
   aryutt[51]=aryutt[6]*aryutt[51];
   aryutt[49]=1./3.*aryutt[51] + 2*aryutt[50] + 2*aryutt[49] - 
   aryutt[23];
   aryutt[50]= - aryutt[6]*aryutt[3];
   aryutt[51]=aryutt[25]*aryutt[55];
   aryutt[51]=16./9.*aryutt[27] - 32./9.*aryutt[30] + aryutt[51];
   aryutt[51]=MMZ*aryutt[51];
   aryutt[38]=4./3.*aryutt[51] + 4*aryutt[50] + 64./27.*aryutt[14] + 47.
   /6.*aryutt[19] - 64./9.*aryutt[33] - 128./27.*aryutt[17] + 92./27.*
   aryutt[35] + 19./6.*aryutt[34] + aryutt[38] + 2*aryutt[24] + 211./54.
    - 2*aryutt[31];
   aryutt[38]=MMZ*aryutt[38];
   aryutt[50]=68*aryutt[17] - 163./2. - 68*aryutt[35];
   aryutt[50]= - 95./6.*aryutt[14] + 1./3.*aryutt[50] + 59./2.*
   aryutt[33];
   aryutt[51]=2*aryutt[30] - aryutt[27];
   aryutt[51]=MMZ*aryutt[51];
   aryutt[50]=34./27.*aryutt[51] + 1./9.*aryutt[50] + 2*aryutt[53];
   aryutt[50]=MMZ*aryutt[50];
   aryutt[51]=17 - 7*aryutt[14];
   aryutt[51]=7./54.*aryutt[51] + aryutt[56];
   aryutt[51]=aryutt[6]*aryutt[51];
   aryutt[53]=77./2. - 34./3.*aryutt[14];
   aryutt[53]=aryutt[15]*aryutt[53];
   aryutt[53]=aryutt[53] - 23*aryutt[23] - 95./6.*aryutt[10] + 68*
   aryutt[21];
   aryutt[50]=aryutt[50] + 1./9.*aryutt[53] + aryutt[51];
   aryutt[50]=aryutt[46]*aryutt[50];
   aryutt[37]=aryutt[37] + aryutt[48] + aryutt[41] + aryutt[50] + 32./9.
   *aryutt[49] + aryutt[38];
   aryutt[37]=aryutt[8]*aryutt[37];
   aryutt[38]=11*aryutt[14] - 7*aryutt[33] - 4*aryutt[17] - 22 - 7*
   aryutt[35];
   aryutt[41]= - 2./3.*aryutt[23] - aryutt[10] + 7./3.*aryutt[21];
   aryutt[41]=aryutt[36]*aryutt[41];
   aryutt[42]=5 + aryutt[42];
   aryutt[42]=aryutt[15]*aryutt[36]*aryutt[42];
   aryutt[48]=aryutt[28] - aryutt[27];
   aryutt[48]=MMt*aryutt[48];
   aryutt[38]=4*aryutt[54] + 8./9.*aryutt[48] + 2./3.*aryutt[42] + 1./9.
   *aryutt[38] + 2*aryutt[41];
   aryutt[41]= - 4*aryutt[14] - 13 + 4*aryutt[17];
   aryutt[41]=aryutt[36]*aryutt[41];
   aryutt[41]=16./9.*aryutt[41] + 2./3.*aryutt[3] - 128./27.*aryutt[30]
    - aryutt[25];
   aryutt[41]=MMZ*aryutt[41];
   aryutt[38]=16./3.*aryutt[38] + aryutt[41];
   aryutt[41]=2./3.*aryutt[14] + 65./3.*aryutt[33] + 14./3.*aryutt[17]
    - 47./8. - 28./3.*aryutt[35];
   aryutt[42]=2./3.*aryutt[23] + aryutt[10] + aryutt[47];
   aryutt[42]=aryutt[36]*aryutt[42];
   aryutt[41]=1./3.*aryutt[41] + 25*aryutt[42];
   aryutt[42]= - 1./3.*aryutt[3] + 14./27.*aryutt[30] + aryutt[27];
   aryutt[42]=2*aryutt[42] + 25./9.*aryutt[52];
   aryutt[42]=MMZ*aryutt[42];
   aryutt[44]= - aryutt[3] + 25./9.*aryutt[44];
   aryutt[44]=aryutt[15]*aryutt[44];
   aryutt[47]=MMt*aryutt[27];
   aryutt[41]=aryutt[42] + 50./3.*aryutt[45] + 28./27.*aryutt[47] + 1./
   3.*aryutt[41] + aryutt[44];
   aryutt[41]=aryutt[46]*aryutt[41];
   aryutt[37]=aryutt[37] + aryutt[39] + aryutt[43] + 2*aryutt[38] + 
   aryutt[41];
   aryutt[37]=aryutt[7]*aryutt[37];
   aryutt[38]=aryutt[3]*aryutt[1];
   aryutt[39]=MMt*aryutt[38];
   aryutt[41]= - aryutt[6]*aryutt[3]*aryutt[1];
   aryutt[39]=aryutt[39] + 9*aryutt[41];
   aryutt[39]=aryutt[6]*aryutt[39];
   aryutt[38]=aryutt[38]*pow(MMt,2);
   aryutt[38]=8*aryutt[38] + aryutt[39];
   aryutt[39]=aryutt[46]*aryutt[38];
   aryutt[38]=aryutt[40]*aryutt[38];
   aryutt[38]=aryutt[39] + aryutt[38];
   aryutt[38]=aryutt[4]*aryutt[38];

      yuttret = aryutt[37] + 4*aryutt[38];
      return yuttret;
}
