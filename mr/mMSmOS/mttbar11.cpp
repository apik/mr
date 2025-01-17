#include <tt.hpp>
namespace mr
{
  double tt<OS>::x11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> armttbar[82], mttbarret;

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
    armttbar[37]=pow(armttbar[2],2);
    armttbar[38]=pow(armttbar[5],2);
    armttbar[39]=armttbar[37] + armttbar[38];
    armttbar[40]=armttbar[39]*armttbar[3];
    armttbar[41]= - 35./9.*armttbar[37] + 104./9. - armttbar[38];
    armttbar[41]=armttbar[4]*armttbar[41];
    armttbar[42]=1./2.*armttbar[38];
    armttbar[43]=17./18.*armttbar[37] + armttbar[42] - 16./9.;
    armttbar[44]=armttbar[8]*armttbar[43];
    armttbar[41]= - 5*armttbar[44] + 3*armttbar[40] + 4./3.*armttbar[41]
      ;
    armttbar[41]=armttbar[8]*armttbar[41];
    armttbar[44]=armttbar[39]*armttbar[4];
    armttbar[45]=25./9.*armttbar[37] + armttbar[38] - 64./9.;
    armttbar[45]=armttbar[45]*armttbar[36];
    armttbar[46]=4*armttbar[45] + armttbar[44];
    armttbar[46]=armttbar[4]*armttbar[46];
    armttbar[41]=armttbar[46] + armttbar[41];
    armttbar[41]=armttbar[6]*armttbar[41];
    armttbar[46]=17./9.*armttbar[37];
    armttbar[47]=armttbar[46] + armttbar[38] - 32./9.;
    armttbar[48]=1./3.*armttbar[8];
    armttbar[49]=armttbar[47]*armttbar[48];
    armttbar[49]=armttbar[49] - armttbar[45];
    armttbar[50]= - armttbar[15]*armttbar[4];
    armttbar[50]=4*armttbar[50] - 2*armttbar[14];
    armttbar[49]=armttbar[49]*armttbar[50];
    armttbar[50]=1./3.*MMZ;
    armttbar[51]=pow(armttbar[8],2);
    armttbar[50]=armttbar[50]*armttbar[51];
    armttbar[52]=armttbar[47]*armttbar[50];
    armttbar[53]=1./2.*armttbar[44];
    armttbar[54]=77./18.*armttbar[37] - 64./9. + 5./2.*armttbar[38];
    armttbar[54]=armttbar[8]*armttbar[54];
    armttbar[41]=armttbar[41] + armttbar[52] + armttbar[54] - 
      armttbar[53] - 5*armttbar[45] - armttbar[40] + armttbar[49];
    armttbar[41]=armttbar[15]*armttbar[41];
    armttbar[49]=armttbar[39]*MMt;
    armttbar[54]=pow(Pi,2);
    armttbar[55]=armttbar[54]*armttbar[49];
    armttbar[56]=armttbar[37] + 1;
    armttbar[56]=armttbar[56]*armttbar[37];
    armttbar[56]=armttbar[56] + armttbar[38];
    armttbar[57]=armttbar[4]*MMt;
    armttbar[58]=armttbar[57]*armttbar[22]*armttbar[56];
    armttbar[55]=armttbar[55] - armttbar[58];
    armttbar[58]=29./4. + 8*armttbar[16];
    armttbar[59]=armttbar[19] + 1;
    armttbar[60]=armttbar[19]*armttbar[59];
    armttbar[58]=5./6.*armttbar[60] + 3*armttbar[35] - 1./3.*
      armttbar[58];
    armttbar[58]=armttbar[39]*armttbar[58];
    armttbar[40]=armttbar[20]*armttbar[40];
    armttbar[40]=13./3.*armttbar[40] + armttbar[58];
    armttbar[40]=MMt*armttbar[40];
    armttbar[58]=armttbar[39]*MMH;
    armttbar[60]=4*armttbar[58] - 16*armttbar[49];
    armttbar[60]=armttbar[29]*armttbar[60];
    armttbar[61]=armttbar[35] - armttbar[16];
    armttbar[61]= - 13 - 4*armttbar[61];
    armttbar[61]=armttbar[61]*armttbar[39];
    armttbar[60]=armttbar[60] + armttbar[61];
    armttbar[60]=MMH*armttbar[60];
    armttbar[61]=2*armttbar[49] - 1./3.*armttbar[58];
    armttbar[61]=MMH*armttbar[61];
    armttbar[62]=pow(MMt,2);
    armttbar[63]=armttbar[39]*armttbar[62];
    armttbar[61]= - 8./3.*armttbar[63] + armttbar[61];
    armttbar[61]=armttbar[26]*armttbar[61];
    armttbar[63]=17./6.*armttbar[22] + 3*armttbar[20];
    armttbar[63]=armttbar[63]*armttbar[39];
    armttbar[64]=armttbar[39]*armttbar[23];
    armttbar[40]=2*armttbar[61] + 1./3.*armttbar[60] + 59./6.*
      armttbar[64] + armttbar[40] + armttbar[63] - 5./6.*armttbar[55];
    armttbar[40]=armttbar[4]*armttbar[40];
    armttbar[55]=52*armttbar[49] - 17*armttbar[58];
    armttbar[55]=armttbar[4]*armttbar[55];
    armttbar[60]=1./2.*armttbar[39];
    armttbar[61]=pow(MMH,2);
    armttbar[63]=armttbar[8]*armttbar[4];
    armttbar[65]=armttbar[61]*armttbar[63]*armttbar[60];
    armttbar[66]=2*armttbar[44];
    armttbar[67]=1./2.*armttbar[58];
    armttbar[68]=armttbar[67]*armttbar[63];
    armttbar[69]=armttbar[66] - armttbar[68];
    armttbar[69]=armttbar[6]*armttbar[69];
    armttbar[66]= - armttbar[13]*armttbar[66];
    armttbar[55]=armttbar[66] + 7*armttbar[69] + armttbar[55] + 
      armttbar[65];
    armttbar[55]=armttbar[12]*armttbar[55];
    armttbar[65]=armttbar[39]*armttbar[57];
    armttbar[66]=armttbar[51]*MMZ;
    armttbar[69]= - 19./2.*armttbar[8] + armttbar[66];
    armttbar[70]=armttbar[38] - 1;
    armttbar[69]=MMZ*armttbar[70]*armttbar[69];
    armttbar[71]=2*armttbar[38];
    armttbar[69]=armttbar[69] + armttbar[71] + 5./2.*armttbar[65];
    armttbar[69]=armttbar[34]*armttbar[69];
    armttbar[68]= - 11*armttbar[44] + armttbar[68];
    armttbar[68]=armttbar[9]*armttbar[68];
    armttbar[72]=armttbar[38]*armttbar[8];
    armttbar[73]= - armttbar[44] - armttbar[72];
    armttbar[74]=armttbar[70]*armttbar[66];
    armttbar[73]=13./2.*armttbar[73] + armttbar[74];
    armttbar[73]=armttbar[11]*armttbar[73];
    armttbar[74]=5*armttbar[38];
    armttbar[75]=armttbar[74]*armttbar[19];
    armttbar[76]= - armttbar[38] + armttbar[75];
    armttbar[76]=armttbar[19]*armttbar[76];
    armttbar[55]=armttbar[73] + armttbar[76] + armttbar[69] + 
      armttbar[55] + armttbar[68];
    armttbar[68]= - armttbar[10]*armttbar[47];
    armttbar[69]= - 34./9.*armttbar[37] + 37./9. + armttbar[38];
    armttbar[69]=armttbar[23]*armttbar[69];
    armttbar[68]=armttbar[69] + armttbar[68];
    armttbar[68]=armttbar[68]*armttbar[48];
    armttbar[69]=armttbar[54]*armttbar[70];
    armttbar[73]=2*armttbar[31];
    armttbar[76]= - 1 + armttbar[73];
    armttbar[76]=armttbar[76]*armttbar[38];
    armttbar[77]= - 163./2. - 68*armttbar[35];
    armttbar[77]=armttbar[77]*armttbar[37];
    armttbar[78]=armttbar[70]*armttbar[19];
    armttbar[68]=armttbar[68] + 5./2.*armttbar[69] - 47./6.*armttbar[78]
      + 1./27.*armttbar[77] + armttbar[76] + 92./27.*armttbar[35] + 211./
      54. - armttbar[73];
    armttbar[68]=armttbar[8]*armttbar[68];
    armttbar[73]=23 - 17*armttbar[37];
    armttbar[69]= - 5*armttbar[69] + 1./9.*armttbar[73] + armttbar[78];
    armttbar[69]=armttbar[8]*armttbar[69];
    armttbar[73]= - armttbar[27]*armttbar[47];
    armttbar[76]=2*armttbar[70];
    armttbar[77]=armttbar[76]*armttbar[25];
    armttbar[73]=armttbar[73] - armttbar[77];
    armttbar[69]=2*armttbar[73] + armttbar[69];
    armttbar[69]=MMZ*armttbar[69]*armttbar[48];
    armttbar[73]=2*armttbar[39];
    armttbar[73]=armttbar[27]*armttbar[73];
    armttbar[79]=25./18.*armttbar[37] - 32./9. + armttbar[42];
    armttbar[79]=armttbar[36]*armttbar[79];
    armttbar[80]= - 1./3.*armttbar[37] + 2./3. - armttbar[38];
    armttbar[81]=2*armttbar[3];
    armttbar[80]=armttbar[80]*armttbar[81];
    armttbar[68]=armttbar[69] + armttbar[68] + armttbar[77] + 
      armttbar[80] + armttbar[73] + 13*armttbar[79];
    armttbar[68]=MMZ*armttbar[68];
    armttbar[69]=armttbar[19] - 1;
    armttbar[69]=armttbar[57]*armttbar[56]*armttbar[69];
    armttbar[60]=armttbar[19]*armttbar[60];
    armttbar[60]=1./6.*armttbar[69] - 1./3.*armttbar[39] + armttbar[60];
    armttbar[60]=armttbar[4]*armttbar[60];
    armttbar[63]=armttbar[39]*armttbar[63];
    armttbar[56]=armttbar[56]*pow(armttbar[4],2);
    armttbar[56]=1./2.*armttbar[56] + armttbar[63];
    armttbar[56]=armttbar[18]*armttbar[56];
    armttbar[69]=11./2.*armttbar[38] + armttbar[75];
    armttbar[69]=armttbar[69]*armttbar[48];
    armttbar[50]= - armttbar[70]*armttbar[50];
    armttbar[70]=armttbar[71]*armttbar[3];
    armttbar[53]=armttbar[70] - armttbar[53];
    armttbar[53]=3*armttbar[53] - 5./3.*armttbar[72];
    armttbar[53]=armttbar[6]*armttbar[8]*armttbar[53];
    armttbar[50]=5./3.*armttbar[56] + armttbar[53] + armttbar[50] + 
      armttbar[69] - armttbar[70] + 5*armttbar[60];
    armttbar[50]=armttbar[18]*armttbar[50];
    armttbar[53]=armttbar[4]*MMH;
    armttbar[56]=armttbar[39]*armttbar[53];
    armttbar[60]=armttbar[56] + armttbar[38] + 41./9.*armttbar[37];
    armttbar[60]=armttbar[8]*armttbar[60];
    armttbar[69]=19./3.*armttbar[37] - 64./3. + armttbar[38];
    armttbar[69]=armttbar[4]*armttbar[69];
    armttbar[60]=armttbar[69] + armttbar[60];
    armttbar[60]=armttbar[8]*armttbar[60];
    armttbar[69]=3*armttbar[38];
    armttbar[70]=25./3.*armttbar[37] + armttbar[69] - 64./3.;
    armttbar[70]=armttbar[70]*armttbar[36];
    armttbar[71]=armttbar[4]*armttbar[70];
    armttbar[72]=MMZ*armttbar[47]*pow(armttbar[8],3);
    armttbar[60]=armttbar[72] - 4*armttbar[71] + armttbar[60];
    armttbar[60]=armttbar[6]*armttbar[60];
    armttbar[46]= - armttbar[78] + armttbar[46] + 22./9. - armttbar[74];
    armttbar[46]=armttbar[46]*armttbar[48];
    armttbar[69]=armttbar[37] - 2 + armttbar[69];
    armttbar[69]=armttbar[3]*armttbar[69];
    armttbar[46]=armttbar[69] + armttbar[46];
    armttbar[69]=MMZ*armttbar[8];
    armttbar[46]=armttbar[46]*armttbar[69];
    armttbar[46]=armttbar[70] - armttbar[46];
    armttbar[71]=4./3.*armttbar[19];
    armttbar[72]=armttbar[71] + 3./2.;
    armttbar[72]=armttbar[38]*armttbar[72];
    armttbar[72]= - 1./3.*armttbar[56] + 119./54.*armttbar[37] + 320./27.
      + armttbar[72];
    armttbar[72]=armttbar[8]*armttbar[72];
    armttbar[71]=23./2. + armttbar[71];
    armttbar[71]=armttbar[71]*armttbar[44];
    armttbar[46]=armttbar[60] + armttbar[72] + armttbar[71] - 2*
      armttbar[46];
    armttbar[46]=armttbar[6]*armttbar[46];
    armttbar[60]=2*armttbar[45];
    armttbar[71]= - 95./18.*armttbar[37] + 64./9. - 7./2.*armttbar[38];
    armttbar[71]=armttbar[71]*armttbar[48];
    armttbar[71]= - armttbar[52] + armttbar[60] + armttbar[71];
    armttbar[71]=MMZ*armttbar[71];
    armttbar[72]=1./9.*armttbar[37];
    armttbar[73]=armttbar[72] + 176./9. + armttbar[74];
    armttbar[73]=1./3.*armttbar[73] - armttbar[65];
    armttbar[42]= - 7./18.*armttbar[37] + 32./9. + armttbar[42];
    armttbar[42]=armttbar[8]*armttbar[42];
    armttbar[43]= - armttbar[43]*armttbar[66];
    armttbar[42]=armttbar[42] + armttbar[43];
    armttbar[42]=armttbar[6]*armttbar[42];
    armttbar[42]=7./3.*armttbar[42] + 2*armttbar[73] + armttbar[71];
    armttbar[42]=armttbar[14]*armttbar[42];
    armttbar[43]= - 4./3.*armttbar[13] - 13./3.*armttbar[6];
    armttbar[66]=armttbar[4]*armttbar[3];
    armttbar[43]=armttbar[39]*armttbar[66]*armttbar[43];
    armttbar[48]= - armttbar[48]*armttbar[56];
    armttbar[48]=7*armttbar[44] + armttbar[48];
    armttbar[43]=1./2.*armttbar[48] + armttbar[43];
    armttbar[43]=armttbar[13]*armttbar[43];
    armttbar[48]=65./9.*armttbar[37] + armttbar[38] - 224./9.;
    armttbar[48]= - armttbar[65] + 1./3.*armttbar[48];
    armttbar[56]=59./18.*armttbar[37] - 64./9. + 3./2.*armttbar[38];
    armttbar[56]=armttbar[8]*armttbar[56];
    armttbar[52]=armttbar[56] - armttbar[52];
    armttbar[52]=MMZ*armttbar[52];
    armttbar[52]=armttbar[52] + armttbar[48];
    armttbar[52]=armttbar[33]*armttbar[52];
    armttbar[48]=armttbar[4]*armttbar[48];
    armttbar[56]=armttbar[47]*armttbar[8];
    armttbar[48]=4*armttbar[56] - 7*armttbar[45] + armttbar[48];
    armttbar[48]=armttbar[21]*armttbar[48];
    armttbar[65]= - 7./9.*armttbar[37] + armttbar[38] + 64./9.;
    armttbar[47]=armttbar[47]*armttbar[69];
    armttbar[47]=armttbar[47] - armttbar[65];
    armttbar[47]=armttbar[30]*MMZ*armttbar[47];
    armttbar[71]= - armttbar[27]*armttbar[65];
    armttbar[71]=64./9.*armttbar[28] + armttbar[71];
    armttbar[71]=MMt*armttbar[71];
    armttbar[47]=armttbar[71] + armttbar[47];
    armttbar[64]= - 5*armttbar[64] + armttbar[67];
    armttbar[53]=armttbar[64]*armttbar[53];
    armttbar[64]= - 23./3.*armttbar[37] - 32./3. + 23./2.*armttbar[38];
    armttbar[64]=armttbar[23]*armttbar[64];
    armttbar[53]=armttbar[64] + armttbar[53];
    armttbar[64]= - 7./6.*armttbar[10] - armttbar[22];
    armttbar[64]=armttbar[64]*armttbar[38];
    armttbar[37]= - 95./54.*armttbar[37] + 64./27.;
    armttbar[37]=armttbar[10]*armttbar[37];
    armttbar[37]=armttbar[64] + armttbar[37] + 1./3.*armttbar[53];
    armttbar[37]=armttbar[8]*armttbar[37];
    armttbar[49]=13./3.*armttbar[49] - 2*armttbar[58];
    armttbar[49]=armttbar[4]*armttbar[49];
    armttbar[53]=armttbar[61]*armttbar[63];
    armttbar[49]=armttbar[49] + 1./6.*armttbar[53];
    armttbar[49]=armttbar[32]*armttbar[49];
    armttbar[45]= - armttbar[45] + 2./3.*armttbar[56];
    armttbar[45]=MMZ*armttbar[45];
    armttbar[45]= - 1./3.*armttbar[65] + armttbar[45];
    armttbar[45]=armttbar[17]*armttbar[45];
    armttbar[53]= - 119./24. + 2*armttbar[35];
    armttbar[38]=armttbar[53]*armttbar[38];
    armttbar[53]= - 47./8. - 28./3.*armttbar[35];
    armttbar[53]=armttbar[53]*armttbar[72];
    armttbar[56]=armttbar[10]*armttbar[70];
    armttbar[58]=armttbar[23]*armttbar[60];
    armttbar[44]=armttbar[25]*armttbar[62]*armttbar[44];
    armttbar[54]= - armttbar[59] + 5*armttbar[54];
    armttbar[54]=armttbar[8]*armttbar[54];
    armttbar[54]=4*armttbar[25] + armttbar[54];
    armttbar[54]=armttbar[8]*armttbar[54];
    armttbar[51]= - armttbar[34]*armttbar[51];
    armttbar[51]=armttbar[54] + armttbar[51];
    armttbar[54]=CW*MMZ;
    armttbar[51]=armttbar[51]*pow(armttbar[54],2);
    armttbar[54]= - armttbar[24]*armttbar[69]*armttbar[76];
    armttbar[59]= - 22 - 7*armttbar[35];
    armttbar[37]=armttbar[41] + 1./3.*armttbar[51] + armttbar[50] + 
      armttbar[48] + armttbar[52] + armttbar[54] + 2*armttbar[45] + 
      armttbar[43] + armttbar[42] + armttbar[49] + armttbar[46] + 
      armttbar[68] + armttbar[37] - 2./3.*armttbar[44] + armttbar[58] + 
      armttbar[56] + armttbar[53] + 32./27.*armttbar[59] + armttbar[38] + 
      armttbar[40] + 4./3.*armttbar[47] + 1./3.*armttbar[55];
    armttbar[37]=armttbar[7]*armttbar[37];
    armttbar[38]=armttbar[62]*armttbar[66];
    armttbar[40]=armttbar[3]*armttbar[57];
    armttbar[41]=armttbar[6]*armttbar[66];
    armttbar[40]=armttbar[40] - 9*armttbar[41];
    armttbar[40]=armttbar[6]*armttbar[40];
    armttbar[38]=8*armttbar[38] + armttbar[40];
    armttbar[38]=armttbar[38]*armttbar[1]*armttbar[39];

    mttbarret = armttbar[37] + 4*armttbar[38];
    return mttbarret.real();
  }
} // namespace mr
