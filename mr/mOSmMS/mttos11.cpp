#include <tt.hpp>
namespace mr
{
  double tt<MS>::x11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> armttos[79], mttosret;

    armttos[1]=double(nH);
    armttos[2]=pow(mmZ,-1);
    armttos[3]=pow(mmH,-1);
    armttos[4]=pow(s,-1);
    armttos[5]=pow(c,-1);
    armttos[6]=Tsil::A(mmt,mu2);
    armttos[7]=double(boson);
    armttos[8]=pow(mmt,-1);
    armttos[9]=Tsil::I2(mmZ,mmt,mmt,mu2);
    armttos[10]=Tsil::I2(mmH,mmt,mmt,mu2);
    armttos[11]=Tsil::I2(0,mmW,mmt,mu2);
    armttos[12]=Tsil::B(mmZ,mmt,mmt,mu2);
    armttos[13]=Tsil::A(mmZ,mu2);
    armttos[14]=Tsil::B(mmH,mmt,mmt,mu2);
    armttos[15]=Tsil::A(mmH,mu2);
    armttos[16]=Tsil::Beps(mmZ,mmt,mmt,mu2);
    armttos[17]=Tsil::Beps(mmH,mmt,mmt,mu2);
    armttos[18]=Tsil::A(mmW,mu2);
    armttos[19]=std::real(Tsil::B(0,mmW,mmt,mu2));
    armttos[20]=Tsil::Aeps(mmW,mu2);
    armttos[21]=Tsil::Aeps(mmZ,mu2);
    armttos[22]=Tsil::Aeps(mmH,mu2);
    armttos[23]=Tsil::Aeps(mmt,mu2);
    armttos[24]=std::real(Tsil::Beps(0,mmW,mmt,mu2));
    armttos[25]=protWt000->M(0);
    armttos[26]=prot0ttHt->M(0);
    armttos[27]=prot0ttZt->M(0);
    armttos[28]=prot0tt0t->M(0);
    armttos[29]=prottH0H->Vxzuv(0);
    armttos[30]=prottZ0Z->Vxzuv(0);
    armttos[31]=protWt000->Uzxyv(0);
    armttos[32]=prot0ttHt->Tuxv(0);
    armttos[33]=prot0ttZt->Tuxv(0);
    armttos[34]=protWt000->Svyz(0);
    armttos[35]=protWt000->Suxv(0);
    armttos[36]=1/(4*mmt - mmZ);
    armttos[37]=pow(armttos[5],2);
    armttos[38]=pow(armttos[4],2);
    armttos[39]=65./9.*armttos[37] + armttos[38] - 224./9.;
    armttos[40]=armttos[33]*armttos[39];
    armttos[41]=armttos[38]*armttos[19];
    armttos[42]=31./4.*armttos[38] - 5*armttos[41];
    armttos[42]=armttos[19]*armttos[42];
    armttos[40]=armttos[42] - armttos[40] - 41./24.*armttos[37] + 460./9.
      + 329./16.*armttos[38];
    armttos[42]=mmH*armttos[26];
    armttos[43]=23./2. + 2*armttos[42];
    armttos[43]=2*armttos[32] + 1./3.*armttos[43];
    armttos[43]=armttos[43]*mmH;
    armttos[44]=13./6.*armttos[11];
    armttos[45]=pow(mmH,2);
    armttos[46]=armttos[45]*armttos[29];
    armttos[43]=armttos[43] - 4./3.*armttos[46] - 3*armttos[22] + 
      armttos[44] + 2./3.*armttos[35] - 17./6.*armttos[20] + 11./3.*
      armttos[10];
    armttos[46]=armttos[37] + armttos[38];
    armttos[43]=armttos[46]*armttos[43];
    armttos[47]=mmH*armttos[46];
    armttos[48]=armttos[47]*armttos[14];
    armttos[49]=armttos[17]*armttos[47];
    armttos[43]=14./3.*armttos[48] - 4./3.*armttos[49] + armttos[43];
    armttos[43]=armttos[2]*armttos[43];
    armttos[49]=armttos[46]*armttos[2];
    armttos[50]=armttos[49]*armttos[19];
    armttos[51]=armttos[37] + 1;
    armttos[51]=armttos[51]*armttos[37];
    armttos[51]=armttos[51] + armttos[38];
    armttos[52]=pow(armttos[2],2);
    armttos[53]=armttos[51]*armttos[52];
    armttos[54]=5./6.*armttos[18];
    armttos[55]= - armttos[53]*armttos[54];
    armttos[56]=armttos[38]*armttos[3];
    armttos[55]=armttos[55] + 10*armttos[56] + armttos[49] - 5./2.*
      armttos[50];
    armttos[55]=armttos[18]*armttos[55];
    armttos[57]=25./9.*armttos[36];
    armttos[57]=armttos[57]*armttos[37];
    armttos[58]=armttos[38]*armttos[36];
    armttos[57]= - 64./9.*armttos[36] + armttos[57] + armttos[58];
    armttos[59]= - 2*armttos[57] - 59./6.*armttos[49];
    armttos[59]=armttos[23]*armttos[59];
    armttos[60]=1./3.*armttos[2];
    armttos[39]= - armttos[39]*armttos[60];
    armttos[39]=7*armttos[57] + armttos[39];
    armttos[39]=armttos[21]*armttos[39];
    armttos[61]=armttos[46]*armttos[14];
    armttos[62]= - 5./2.*armttos[46] + 2./3.*armttos[61];
    armttos[62]=armttos[2]*armttos[62];
    armttos[63]=armttos[46]*armttos[3]*armttos[2];
    armttos[64]=armttos[63]*armttos[15];
    armttos[62]=armttos[62] + 4./3.*armttos[64];
    armttos[62]=armttos[15]*armttos[62];
    armttos[65]= - armttos[38] - 64./9. + 7./9.*armttos[37];
    armttos[66]=2./3.*armttos[65];
    armttos[67]= - armttos[16]*armttos[66];
    armttos[68]=5./2.*armttos[38];
    armttos[69]= - 13./18.*armttos[37] - 64./9. - armttos[68];
    armttos[69]=armttos[12]*armttos[69];
    armttos[70]=armttos[37]*armttos[36];
    armttos[71]= - 64./3.*armttos[36] + 25./3.*armttos[70] + 3*
      armttos[58];
    armttos[72]= - armttos[9]*armttos[71];
    armttos[39]=armttos[55] + armttos[62] + armttos[72] + armttos[69] + 
      armttos[39] + armttos[67] + armttos[59] + armttos[43] + 1./3.*
      armttos[40];
    armttos[39]=armttos[7]*armttos[39];
    armttos[40]= - armttos[31] + armttos[24];
    armttos[43]=armttos[38] - 1;
    armttos[55]=2*armttos[43];
    armttos[40]=armttos[55]*armttos[40];
    armttos[59]= - 59./18.*armttos[37] + 64./9. - 3./2.*armttos[38];
    armttos[59]=armttos[33]*armttos[59];
    armttos[62]=17./9.*armttos[37] + armttos[38] - 32./9.;
    armttos[67]=armttos[16]*armttos[62];
    armttos[69]=1./2.*armttos[38];
    armttos[72]= - 41./18.*armttos[37] + 64./9. - armttos[69];
    armttos[72]=armttos[12]*armttos[72];
    armttos[73]=pow(Pi,2);
    armttos[74]=armttos[43]*armttos[73];
    armttos[75]=armttos[43]*armttos[19];
    armttos[76]=161./27. - 19*armttos[38];
    armttos[40]= - 5./2.*armttos[74] + 1./3.*armttos[72] + 41./12.*
      armttos[75] - 4./3.*armttos[67] + armttos[59] + 1./8.*armttos[76] + 
      56./27.*armttos[37] + armttos[40];
    armttos[40]=armttos[7]*armttos[40];
    armttos[59]= - 1./2.*armttos[18] - armttos[11] + armttos[35];
    armttos[59]=armttos[43]*armttos[59];
    armttos[67]=34./9.*armttos[37];
    armttos[72]=armttos[67] - 37./9. - armttos[38];
    armttos[72]=armttos[23]*armttos[72];
    armttos[74]=armttos[9]*armttos[62];
    armttos[59]=armttos[74] + armttos[72] + armttos[59];
    armttos[59]=armttos[7]*armttos[59];
    armttos[72]=pow(c,2);
    armttos[43]=armttos[72] - armttos[43];
    armttos[72]= - 5*armttos[73] + 1./2.*armttos[19];
    armttos[72]=armttos[43]*armttos[72];
    armttos[74]=armttos[12] + armttos[33] + 1;
    armttos[74]=armttos[62]*armttos[74];
    armttos[72]=armttos[74] + armttos[72];
    armttos[74]=mmZ*armttos[7];
    armttos[72]=armttos[72]*armttos[74];
    armttos[59]=armttos[59] + armttos[72];
    armttos[72]=1./3.*armttos[8];
    armttos[59]=armttos[59]*armttos[72];
    armttos[43]= - armttos[25]*armttos[43];
    armttos[76]= - armttos[30]*armttos[62];
    armttos[43]=armttos[43] + armttos[76];
    armttos[43]=armttos[43]*armttos[74];
    armttos[40]=armttos[59] + 4./3.*armttos[43] + armttos[40];
    armttos[40]=mmZ*armttos[40];
    armttos[43]= - armttos[14] - armttos[32] - 1;
    armttos[43]=armttos[2]*armttos[43]*armttos[46]*armttos[45];
    armttos[45]=armttos[47]*armttos[2];
    armttos[59]=armttos[10]*armttos[45];
    armttos[43]=armttos[43] - armttos[59];
    armttos[59]=23./2.*armttos[38];
    armttos[76]=5*armttos[45] + 23./3.*armttos[37] + 32./3. - 
      armttos[59];
    armttos[76]=armttos[23]*armttos[76];
    armttos[77]=95./18.*armttos[37] - 64./9. + 7./2.*armttos[38];
    armttos[77]=armttos[9]*armttos[77];
    armttos[43]=armttos[77] + armttos[76] + 1./2.*armttos[43];
    armttos[44]=armttos[44] - 19./6.*armttos[35] + armttos[20];
    armttos[44]=armttos[38]*armttos[44];
    armttos[76]=armttos[21]*armttos[62];
    armttos[77]=armttos[15]*armttos[45];
    armttos[78]=armttos[18]*armttos[49];
    armttos[78]=armttos[41] + armttos[78];
    armttos[78]=1./4.*armttos[38] - 5./3.*armttos[78];
    armttos[78]=armttos[18]*armttos[78];
    armttos[43]=armttos[78] + 1./6.*armttos[77] - 4*armttos[76] + 
      armttos[44] + 1./3.*armttos[43];
    armttos[43]=armttos[7]*armttos[43];
    armttos[40]=armttos[43] + armttos[40];
    armttos[40]=armttos[8]*armttos[40];
    armttos[43]=5./6.*armttos[51];
    armttos[44]=armttos[35] - armttos[20];
    armttos[43]=armttos[43]*armttos[44];
    armttos[44]=armttos[21]*armttos[46];
    armttos[43]=armttos[44] + armttos[43];
    armttos[43]=armttos[52]*armttos[43];
    armttos[44]=armttos[29]*mmH;
    armttos[42]=8./3.*armttos[17] + 103./24. - 13./3.*armttos[32] + 16./
      3.*armttos[44] - 4*armttos[42];
    armttos[42]=armttos[46]*armttos[42];
    armttos[42]= - 46./3.*armttos[61] + armttos[42];
    armttos[42]=armttos[2]*armttos[42];
    armttos[44]=1./2.*armttos[49];
    armttos[52]= - armttos[44] - armttos[50];
    armttos[52]=armttos[19]*armttos[52];
    armttos[76]=armttos[53]*armttos[19];
    armttos[53]= - 1./2.*armttos[53] - armttos[76];
    armttos[53]=armttos[53]*armttos[54];
    armttos[54]=2*armttos[12];
    armttos[77]=armttos[54] + armttos[33];
    armttos[77]=armttos[49]*armttos[77];
    armttos[78]= - armttos[22]*armttos[63];
    armttos[42]=armttos[53] + 13./3.*armttos[78] + 5./6.*armttos[52] - 
      256./27.*armttos[28] + armttos[42] + armttos[43] + armttos[77];
    armttos[42]=armttos[7]*armttos[42];
    armttos[43]=8*armttos[26] + armttos[25];
    armttos[43]=armttos[46]*armttos[43];
    armttos[51]=armttos[2]*armttos[51];
    armttos[43]=2*armttos[43] - 25./16.*armttos[51];
    armttos[43]=armttos[2]*armttos[43];
    armttos[43]=armttos[43] - 5./4.*armttos[76];
    armttos[43]=armttos[7]*armttos[43];
    armttos[51]=armttos[63]*armttos[1];
    armttos[43]= - 24*armttos[51] + 1./3.*armttos[43];
    armttos[43]=mmt*armttos[43];
    armttos[52]=armttos[7]*armttos[2];
    armttos[53]=armttos[46]*armttos[52];
    armttos[63]=armttos[73]*armttos[53];
    armttos[42]=armttos[43] + armttos[42] + 5./6.*armttos[63];
    armttos[42]=mmt*armttos[42];
    armttos[43]= - 47*armttos[46] + 4*armttos[61];
    armttos[43]=armttos[43]*armttos[60];
    armttos[61]= - armttos[71]*armttos[54];
    armttos[43]=13./3.*armttos[64] + armttos[61] + 2./3.*armttos[50] + 
      10*armttos[57] + armttos[43];
    armttos[43]=armttos[7]*armttos[43];
    armttos[47]=1./2.*armttos[48] + armttos[47];
    armttos[47]=armttos[2]*armttos[47];
    armttos[48]=253./18.*armttos[37] - 416./9. + armttos[68];
    armttos[48]=armttos[12]*armttos[48];
    armttos[41]=armttos[48] + 2*armttos[41] + 7*armttos[47] - 7./18.*
      armttos[37] + 272./9. - 55./2.*armttos[38];
    armttos[47]=armttos[15]*armttos[49];
    armttos[48]=7./2.*armttos[49] + 6*armttos[56];
    armttos[48]=armttos[18]*armttos[48];
    armttos[41]=armttos[48] + 1./3.*armttos[41] + 3*armttos[47];
    armttos[41]=armttos[7]*armttos[41];
    armttos[47]=5*armttos[38];
    armttos[48]=4*armttos[75] + armttos[67] - 91./9. + armttos[47];
    armttos[50]=17./18.*armttos[37] + armttos[69] - 16./9.;
    armttos[56]=armttos[12]*armttos[50];
    armttos[48]=2*armttos[48] + 13*armttos[56];
    armttos[48]=armttos[48]*armttos[74];
    armttos[56]=armttos[7]*armttos[18]*armttos[38];
    armttos[48]=11*armttos[56] + armttos[48];
    armttos[48]=armttos[48]*armttos[72];
    armttos[56]=armttos[37] - 2 + 3*armttos[38];
    armttos[56]=armttos[3]*armttos[56]*armttos[74];
    armttos[41]=armttos[48] + armttos[41] + 2*armttos[56];
    armttos[41]=armttos[8]*armttos[41];
    armttos[48]=25*armttos[70] - 64*armttos[36] + 9*armttos[58];
    armttos[48]=armttos[48]*armttos[52];
    armttos[48]= - 9*armttos[51] + armttos[48];
    armttos[56]= - 29./9.*armttos[37] + 64./9. - armttos[38];
    armttos[56]=2*armttos[56] - armttos[45];
    armttos[56]=armttos[7]*armttos[56];
    armttos[61]=armttos[74]*armttos[8];
    armttos[63]=armttos[61]*armttos[62];
    armttos[56]=armttos[56] - armttos[63];
    armttos[56]=armttos[8]*armttos[56];
    armttos[47]= - 21*armttos[37] + 64 - armttos[47];
    armttos[47]=armttos[47]*armttos[52];
    armttos[47]=armttos[47] + armttos[56];
    armttos[47]=armttos[8]*armttos[47];
    armttos[47]=4*armttos[48] + armttos[47];
    armttos[47]=armttos[6]*armttos[47];
    armttos[48]=mmt*armttos[51];
    armttos[41]=armttos[47] + armttos[41] + armttos[43] - 28*armttos[48]
      ;
    armttos[41]=armttos[6]*armttos[41];
    armttos[43]=armttos[46]*armttos[3];
    armttos[47]= - armttos[57]*armttos[54];
    armttos[44]=5*armttos[43] + armttos[47] + 2*armttos[71] + 
      armttos[44];
    armttos[44]=armttos[7]*armttos[44];
    armttos[47]= - 16*armttos[57] - armttos[49];
    armttos[47]=armttos[47]*armttos[52];
    armttos[48]=365./9.*armttos[37] - 992./9. + 13*armttos[38];
    armttos[48]=armttos[48]*armttos[60];
    armttos[43]=armttos[48] + 3*armttos[43];
    armttos[43]=armttos[7]*armttos[43];
    armttos[48]=armttos[8]*armttos[7];
    armttos[49]=armttos[50]*armttos[48];
    armttos[43]=armttos[43] + 7*armttos[49];
    armttos[43]=armttos[8]*armttos[43];
    armttos[43]=armttos[47] + armttos[43];
    armttos[43]=armttos[6]*armttos[43];
    armttos[47]=armttos[62]*armttos[54];
    armttos[47]=armttos[47] - 367./18.*armttos[37] + 320./9. - 
      armttos[59];
    armttos[47]=armttos[7]*armttos[47];
    armttos[47]=armttos[47] - armttos[63];
    armttos[47]=armttos[47]*armttos[72];
    armttos[49]=armttos[62]*armttos[72];
    armttos[49]= - armttos[57] + armttos[49];
    armttos[49]=armttos[13]*armttos[52]*armttos[49];
    armttos[43]=4*armttos[49] + armttos[43] + armttos[44] + armttos[47];
    armttos[43]=armttos[13]*armttos[43];
    armttos[44]=2*armttos[16] - 7;
    armttos[44]=armttos[57]*armttos[44];
    armttos[47]=armttos[30]*armttos[65];
    armttos[49]= - 25./6.*armttos[70] + 32./3.*armttos[36] - 3./2.*
      armttos[58];
    armttos[49]=armttos[12]*armttos[49];
    armttos[50]=1./3.*armttos[37] - 2./3. + armttos[38];
    armttos[50]=armttos[3]*armttos[50];
    armttos[51]= - armttos[25]*armttos[55];
    armttos[44]=10*armttos[50] + armttos[49] - 4./3.*armttos[47] + 
      armttos[51] + armttos[44];
    armttos[44]=armttos[44]*armttos[74];
    armttos[46]= - armttos[46]*armttos[74];
    armttos[47]= - mmt*armttos[7]*armttos[66];
    armttos[48]=armttos[62]*pow(mmZ,2)*armttos[48];
    armttos[46]=1./3.*armttos[48] + armttos[46] + armttos[47];
    armttos[46]=armttos[27]*armttos[46];
    armttos[38]= - 2./3.*armttos[45] - 14./27.*armttos[37] - 112./27. + 
      armttos[38];
    armttos[38]=armttos[7]*armttos[38];
    armttos[37]=23 - 17*armttos[37];
    armttos[37]=armttos[37]*armttos[61];
    armttos[37]=armttos[38] + 2./27.*armttos[37];
    armttos[37]=armttos[8]*armttos[37];
    armttos[37]=3*armttos[53] + 2*armttos[37];
    armttos[37]=armttos[34]*armttos[37];

    mttosret = armttos[37] + armttos[39] + armttos[40] + armttos[41]
      + armttos[42] + armttos[43] + armttos[44] + 2*armttos[46];
    return mttosret.real();
  }
} // namespace mr
