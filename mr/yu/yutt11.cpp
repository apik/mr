#include <tt.hpp>
namespace mr
{
  double tt<OS>::y11(size_t nL, size_t nH, size_t boson)
  {     
      
      
    std::complex<double> aryutt[73], yuttret;

    aryutt[1]=double(nH);
    aryutt[2]=pow(CW,-1);
    aryutt[3]=pow(MMZ,-1);
    aryutt[4]=pow(SW,-1);
    aryutt[5]=Tsil::I2(0,0,MMt,mu2);
    aryutt[6]=Tsil::A(MMt,mu2);
    aryutt[7]=pow(MMt,-1);
    aryutt[8]=Tsil::Aeps(MMt,mu2);
    aryutt[9]=double(boson);
    aryutt[10]=Tsil::I2(MMH,MMt,MMt,mu2);
    aryutt[11]=Tsil::I2(MMZ,MMt,MMt,mu2);
    aryutt[12]=Tsil::I2(0,MMW,MMt,mu2);
    aryutt[13]=Tsil::B(MMH,MMt,MMt,mu2);
    aryutt[14]=Tsil::A(MMH,mu2);
    aryutt[15]=Tsil::B(MMZ,MMt,MMt,mu2);
    aryutt[16]=Tsil::A(MMZ,mu2);
    aryutt[17]=Tsil::Beps(MMH,MMt,MMt,mu2);
    aryutt[18]=Tsil::Beps(MMZ,MMt,MMt,mu2);
    aryutt[19]=pow(MMH,-1);
    aryutt[20]=Tsil::A(MMW,mu2);
    aryutt[21]=std::real(Tsil::B(0,MMW,MMt,mu2));
    aryutt[22]=Tsil::Aeps(MMH,mu2);
    aryutt[23]=Tsil::Aeps(MMZ,mu2);
    aryutt[24]=Tsil::Aeps(MMW,mu2);
    aryutt[25]=std::real(Tsil::Beps(0,MMW,MMt,mu2));
    aryutt[26]=protWt000->M(0);
    aryutt[27]=prot0ttHt->M(0);
    aryutt[28]=prot0ttZt->M(0);
    aryutt[29]=prot0tt0t->M(0);
    aryutt[30]=prottH0H->Vxzuv(0);
    aryutt[31]=prottZ0Z->Vxzuv(0);
    aryutt[32]=protWt000->Uzxyv(0);
    aryutt[33]=prot0ttHt->Tuxv(0);
    aryutt[34]=prot0ttZt->Tuxv(0);
    aryutt[35]=protWt000->Txuv(0);
    aryutt[36]=protWt000->Tyzv(0);
    aryutt[37]=1/(4*MMt - MMZ);
    aryutt[38]=1/( - MMW + MMH);
    aryutt[39]=pow(aryutt[2],2);
    aryutt[40]=pow(aryutt[4],2);
    aryutt[41]=aryutt[39] + aryutt[40];
    aryutt[42]=aryutt[41]*aryutt[21];
    aryutt[43]=aryutt[41]*aryutt[13];
    aryutt[44]=aryutt[41]*aryutt[3];
    aryutt[45]=aryutt[16]*aryutt[44];
    aryutt[45]=aryutt[45] + 14./3.*aryutt[43] + 23./2.*aryutt[41] + 4./3.
      *aryutt[42];
    aryutt[45]=aryutt[6]*aryutt[45];
    aryutt[46]=aryutt[41]*aryutt[33];
    aryutt[47]=aryutt[46] + aryutt[43] + aryutt[41];
    aryutt[48]=1./2.*aryutt[7];
    aryutt[47]=aryutt[48]*aryutt[47];
    aryutt[49]=aryutt[41]*aryutt[27];
    aryutt[50]=2*aryutt[41];
    aryutt[50]=aryutt[30]*aryutt[50];
    aryutt[50]= - aryutt[49] + aryutt[50];
    aryutt[47]=2*aryutt[50] + aryutt[47];
    aryutt[47]=MMH*aryutt[47];
    aryutt[47]=aryutt[47] - 53./4.*aryutt[41] - 17*aryutt[43];
    aryutt[50]= - 7*aryutt[43] - 1./2.*aryutt[41];
    aryutt[50]=aryutt[6]*aryutt[50];
    aryutt[51]=aryutt[41]*aryutt[10];
    aryutt[50]=aryutt[50] + aryutt[51];
    aryutt[52]=pow(aryutt[6],2);
    aryutt[53]=aryutt[52]*aryutt[7];
    aryutt[54]=aryutt[53] - 1./6.*aryutt[14] - 5./3.*aryutt[8];
    aryutt[54]=aryutt[41]*aryutt[54];
    aryutt[50]=1./6.*aryutt[50] + aryutt[54];
    aryutt[50]=aryutt[7]*aryutt[50];
    aryutt[54]=4./3.*aryutt[41];
    aryutt[55]=aryutt[17]*aryutt[54];
    aryutt[46]=aryutt[55] - 2*aryutt[46] + aryutt[50] + 1./3.*aryutt[47]
      ;
    aryutt[46]=MMH*aryutt[46];
    aryutt[47]=aryutt[7]*aryutt[6];
    aryutt[50]=3./2.*aryutt[47];
    aryutt[55]= - aryutt[50] + 4;
    aryutt[55]=aryutt[41]*aryutt[55];
    aryutt[56]=aryutt[19]*aryutt[41];
    aryutt[57]= - 4./3.*aryutt[14] - 13./3.*aryutt[6];
    aryutt[57]=aryutt[56]*aryutt[57];
    aryutt[43]= - 2./3.*aryutt[43] + aryutt[57] + aryutt[55];
    aryutt[43]=aryutt[14]*aryutt[43];
    aryutt[55]=1./2.*aryutt[40];
    aryutt[57]= - 1 + aryutt[55];
    aryutt[57]=aryutt[57]*aryutt[40];
    aryutt[58]=3*aryutt[40];
    aryutt[59]=1./3. - aryutt[58];
    aryutt[59]=aryutt[59]*aryutt[55];
    aryutt[59]= - 199./54.*aryutt[39] + 416./27. + aryutt[59];
    aryutt[59]=aryutt[59]*aryutt[47];
    aryutt[60]=17./9.*aryutt[39];
    aryutt[61]=aryutt[60] + aryutt[40] - 32./9.;
    aryutt[62]=aryutt[61]*aryutt[16];
    aryutt[63]=aryutt[7]*aryutt[62];
    aryutt[57]= - 4./3.*aryutt[63] + aryutt[59] + aryutt[57] - 
      aryutt[39];
    aryutt[57]=aryutt[16]*aryutt[57];
    aryutt[59]=13./6.*aryutt[12];
    aryutt[63]= - aryutt[59] + 59./6.*aryutt[8] + 3*aryutt[22];
    aryutt[63]=aryutt[41]*aryutt[63];
    aryutt[64]=19./3.*aryutt[39] - 64./3. + aryutt[40];
    aryutt[64]=aryutt[64]*aryutt[53];
    aryutt[65]=aryutt[40] - 224./9. + 65./9.*aryutt[39];
    aryutt[65]=1./3.*aryutt[65];
    aryutt[66]=aryutt[23]*aryutt[65];
    aryutt[43]=aryutt[46] + aryutt[57] + aryutt[43] + aryutt[66] + 
      aryutt[64] - 11./3.*aryutt[51] + aryutt[45] + aryutt[63];
    aryutt[43]=aryutt[3]*aryutt[43];
    aryutt[45]=aryutt[61]*aryutt[15];
    aryutt[46]=aryutt[55] - 16./9. + 17./18.*aryutt[39];
    aryutt[51]=aryutt[46]*aryutt[47];
    aryutt[51]= - 5*aryutt[51] - 2./3.*aryutt[45] + 77./18.*aryutt[39]
      - 64./9. + 5./2.*aryutt[40];
    aryutt[51]=aryutt[16]*aryutt[51];
    aryutt[57]=4./3.*aryutt[21] + 9./4.;
    aryutt[57]=aryutt[40]*aryutt[57];
    aryutt[63]= - 7./18.*aryutt[39] + 32./9. + aryutt[55];
    aryutt[63]=aryutt[15]*aryutt[63];
    aryutt[57]=7./3.*aryutt[63] + 265./108.*aryutt[39] + 320./27. + 
      aryutt[57];
    aryutt[57]=aryutt[6]*aryutt[57];
    aryutt[63]=95./18.*aryutt[39] - 64./9. + 7./2.*aryutt[40];
    aryutt[63]=1./3.*aryutt[63];
    aryutt[64]= - aryutt[11]*aryutt[63];
    aryutt[66]=aryutt[40] + 41./9.*aryutt[39];
    aryutt[66]=aryutt[66]*aryutt[53];
    aryutt[67]=aryutt[23]*aryutt[61];
    aryutt[68]= - 23./3.*aryutt[39] - 32./3. + 23./2.*aryutt[40];
    aryutt[69]=1./3.*aryutt[8];
    aryutt[68]=aryutt[68]*aryutt[69];
    aryutt[59]= - aryutt[40]*aryutt[59];
    aryutt[51]=aryutt[51] + aryutt[59] + aryutt[68] + 4*aryutt[67] + 
      aryutt[66] + aryutt[64] + aryutt[57];
    aryutt[51]=aryutt[7]*aryutt[51];
    aryutt[57]=25./9.*aryutt[39] + aryutt[40] - 64./9.;
    aryutt[59]=2*aryutt[15];
    aryutt[64]=aryutt[57]*aryutt[59];
    aryutt[66]= - 5*aryutt[57] + aryutt[64];
    aryutt[66]=aryutt[16]*aryutt[66];
    aryutt[67]=25./3.*aryutt[39] + aryutt[58] - 64./3.;
    aryutt[52]= - aryutt[67]*aryutt[52];
    aryutt[68]=aryutt[6] + aryutt[16];
    aryutt[68]=aryutt[16]*aryutt[57]*aryutt[68];
    aryutt[52]=aryutt[52] + aryutt[68];
    aryutt[52]=aryutt[3]*aryutt[52];
    aryutt[68]=2*aryutt[8];
    aryutt[70]=aryutt[68] - 7*aryutt[23];
    aryutt[70]=aryutt[57]*aryutt[70];
    aryutt[71]=aryutt[11] - 2*aryutt[6];
    aryutt[67]=aryutt[67]*aryutt[71];
    aryutt[52]=4*aryutt[52] + aryutt[66] + aryutt[67] + aryutt[70];
    aryutt[52]=aryutt[37]*aryutt[52];
    aryutt[66]=2*aryutt[35] - 125./8.;
    aryutt[66]=aryutt[40]*aryutt[66];
    aryutt[67]=5*aryutt[40];
    aryutt[70]=aryutt[67]*aryutt[21];
    aryutt[71]= - aryutt[40] + aryutt[70];
    aryutt[71]=aryutt[21]*aryutt[71];
    aryutt[72]=1./9.*aryutt[39] + 176./9. + aryutt[67];
    aryutt[72]=aryutt[72]*aryutt[59];
    aryutt[66]=aryutt[72] + aryutt[71] - 53./24.*aryutt[39] - 704./9. + 
      aryutt[66];
    aryutt[71]=MMH*aryutt[44];
    aryutt[71]= - 2./3.*aryutt[71] - 14./27.*aryutt[39] - 112./27. + 
      aryutt[40];
    aryutt[71]=aryutt[36]*aryutt[71];
    aryutt[65]=aryutt[34]*aryutt[65];
    aryutt[58]=aryutt[47]*aryutt[58];
    aryutt[58]= - aryutt[40] + aryutt[58];
    aryutt[58]=aryutt[14]*aryutt[38]*aryutt[58];
    aryutt[72]= - aryutt[7]*aryutt[40];
    aryutt[44]=aryutt[72] + 17./6.*aryutt[44];
    aryutt[44]=aryutt[24]*aryutt[44];
    aryutt[43]=aryutt[52] + aryutt[44] + 2*aryutt[71] + aryutt[43] + 1./
      2.*aryutt[58] + aryutt[65] + 1./3.*aryutt[66] + aryutt[51];
    aryutt[43]=aryutt[9]*aryutt[43];
    aryutt[44]=aryutt[40] - 1;
    aryutt[51]=aryutt[44]*aryutt[21];
    aryutt[52]= - aryutt[51] + aryutt[60] + 22./9. - aryutt[67];
    aryutt[46]=aryutt[15]*aryutt[46];
    aryutt[46]= - 7*aryutt[46] + 2*aryutt[52];
    aryutt[46]=aryutt[6]*aryutt[46];
    aryutt[52]= - aryutt[11]*aryutt[61];
    aryutt[46]=aryutt[52] + aryutt[46];
    aryutt[52]=aryutt[61]*aryutt[53];
    aryutt[46]=1./3.*aryutt[46] + aryutt[52];
    aryutt[46]=aryutt[7]*aryutt[46];
    aryutt[52]= - 2*aryutt[25] - 19./6.*aryutt[35];
    aryutt[52]=aryutt[44]*aryutt[52];
    aryutt[58]=2*aryutt[32];
    aryutt[60]= - 1 + aryutt[58];
    aryutt[60]=aryutt[60]*aryutt[40];
    aryutt[63]= - aryutt[15]*aryutt[63];
    aryutt[65]=59./18.*aryutt[39] - 64./9. + 3./2.*aryutt[40];
    aryutt[65]=aryutt[34]*aryutt[65];
    aryutt[66]= - 23 + 17*aryutt[39];
    aryutt[67]=aryutt[36]*aryutt[66];
    aryutt[46]= - 4./27.*aryutt[67] + aryutt[65] + aryutt[46] + 
      aryutt[63] - 47./6.*aryutt[51] - 163./54.*aryutt[39] + aryutt[60] + 
      211./54. - aryutt[58] + aryutt[52];
    aryutt[46]=aryutt[7]*aryutt[46];
    aryutt[52]=aryutt[12]*aryutt[44];
    aryutt[52]=aryutt[62] + aryutt[52];
    aryutt[58]= - 34./9.*aryutt[39] + 37./9. + aryutt[40];
    aryutt[58]=aryutt[58]*aryutt[69];
    aryutt[52]=aryutt[58] + 1./3.*aryutt[52];
    aryutt[58]=pow(aryutt[7],2);
    aryutt[52]=aryutt[58]*aryutt[52];
    aryutt[60]=aryutt[28]*aryutt[41];
    aryutt[62]= - 7./9.*aryutt[39] + aryutt[40] + 64./9.;
    aryutt[63]=aryutt[31]*aryutt[62];
    aryutt[65]=aryutt[26]*aryutt[44];
    aryutt[60]= - 2./3.*aryutt[63] + aryutt[65] + aryutt[60];
    aryutt[63]=25./18.*aryutt[39] - 32./9. + aryutt[55];
    aryutt[63]=13*aryutt[63] + aryutt[64];
    aryutt[63]=aryutt[37]*aryutt[63];
    aryutt[46]=aryutt[63] + 2*aryutt[60] + aryutt[52] + aryutt[46];
    aryutt[46]=aryutt[9]*aryutt[46];
    aryutt[52]=pow(CW,2);
    aryutt[60]= - 1 - aryutt[21];
    aryutt[60]=aryutt[60]*aryutt[52];
    aryutt[52]=aryutt[44] - aryutt[52];
    aryutt[63]=aryutt[35]*aryutt[52];
    aryutt[45]=aryutt[63] - aryutt[45] + aryutt[60] - 1./9.*aryutt[66]
      + aryutt[51];
    aryutt[45]=aryutt[7]*aryutt[45];
    aryutt[51]=2*aryutt[31] - aryutt[28];
    aryutt[51]=aryutt[61]*aryutt[51];
    aryutt[60]= - 2*aryutt[52];
    aryutt[60]=aryutt[26]*aryutt[60];
    aryutt[51]=aryutt[60] + aryutt[51];
    aryutt[45]=2*aryutt[51] + aryutt[45];
    aryutt[45]=aryutt[7]*aryutt[45];
    aryutt[51]= - aryutt[34]*aryutt[61]*aryutt[58];
    aryutt[45]=aryutt[45] + aryutt[51];
    aryutt[51]=1./3.*aryutt[9];
    aryutt[60]=aryutt[51]*MMZ;
    aryutt[45]=aryutt[45]*aryutt[60];
    aryutt[45]=aryutt[45] + aryutt[46];
    aryutt[45]=MMZ*aryutt[45];
    aryutt[46]= - aryutt[3]*aryutt[23];
    aryutt[46]=aryutt[46] - 8./3.*aryutt[17] - aryutt[34] + 13./3.*
      aryutt[33] + 52./3.*aryutt[13] + 5./6.*aryutt[35] - aryutt[59];
    aryutt[46]=aryutt[41]*aryutt[46];
    aryutt[59]=aryutt[42] + aryutt[41];
    aryutt[59]=aryutt[21]*aryutt[59];
    aryutt[59]= - 29./2.*aryutt[41] + 5*aryutt[59];
    aryutt[54]= - aryutt[30]*aryutt[54];
    aryutt[54]=aryutt[49] + aryutt[54];
    aryutt[54]=MMH*aryutt[54];
    aryutt[56]=aryutt[22]*aryutt[56];
    aryutt[46]=4*aryutt[54] + 13./3.*aryutt[56] + 1./6.*aryutt[59] + 
      aryutt[46];
    aryutt[46]=aryutt[3]*aryutt[46];
    aryutt[54]= - aryutt[28]*aryutt[62];
    aryutt[54]=64./9.*aryutt[29] + aryutt[54];
    aryutt[56]=aryutt[39] + 1;
    aryutt[56]=aryutt[56]*aryutt[39];
    aryutt[56]=aryutt[56] + aryutt[40];
    aryutt[59]=pow(aryutt[3],2);
    aryutt[59]=5./6.*aryutt[59];
    aryutt[59]=aryutt[59]*aryutt[56];
    aryutt[63]=aryutt[24]*aryutt[59];
    aryutt[46]=aryutt[63] + 4./3.*aryutt[54] + aryutt[46];
    aryutt[46]=aryutt[9]*aryutt[46];
    aryutt[54]= - aryutt[26]*aryutt[41];
    aryutt[49]= - 8*aryutt[49] + aryutt[54];
    aryutt[49]=MMt*aryutt[49];
    aryutt[54]=aryutt[36]*aryutt[41];
    aryutt[49]=2./3.*aryutt[49] + 3*aryutt[54];
    aryutt[54]=aryutt[9]*aryutt[3];
    aryutt[49]=aryutt[54]*aryutt[49];
    aryutt[63]=aryutt[3]*aryutt[1]*aryutt[41];
    aryutt[46]=39./4.*aryutt[63] + aryutt[49] + aryutt[46];
    aryutt[46]=MMt*aryutt[46];
    aryutt[49]= - 8./3. - aryutt[55];
    aryutt[49]=aryutt[49]*aryutt[40];
    aryutt[64]=1 + aryutt[40];
    aryutt[64]=aryutt[64]*aryutt[40];
    aryutt[64]=aryutt[64] + aryutt[39];
    aryutt[50]=aryutt[64]*aryutt[50];
    aryutt[39]=aryutt[50] + 5./2.*aryutt[42] + aryutt[49] - 8./3.*
      aryutt[39];
    aryutt[39]=aryutt[3]*aryutt[39];
    aryutt[42]=aryutt[6]*aryutt[38];
    aryutt[42]= - 5./3.*aryutt[47] - 3./2.*aryutt[42];
    aryutt[42]=aryutt[40]*aryutt[42];
    aryutt[40]=11./2.*aryutt[40] + aryutt[70];
    aryutt[40]=1./3.*aryutt[40] + aryutt[42];
    aryutt[40]=aryutt[7]*aryutt[40];
    aryutt[42]=aryutt[21] - 1;
    aryutt[42]=MMt*aryutt[42]*aryutt[59];
    aryutt[47]=aryutt[38]*aryutt[55];
    aryutt[39]=aryutt[42] + aryutt[39] + aryutt[47] + aryutt[40];
    aryutt[39]=aryutt[9]*aryutt[39];
    aryutt[40]=aryutt[7]*aryutt[41];
    aryutt[42]=aryutt[3]*aryutt[56];
    aryutt[40]=aryutt[40] + 1./2.*aryutt[42];
    aryutt[40]=aryutt[20]*aryutt[40]*aryutt[54];
    aryutt[42]=aryutt[60]*aryutt[58];
    aryutt[47]= - aryutt[44]*aryutt[42];
    aryutt[39]=5./3.*aryutt[40] + aryutt[47] + aryutt[39];
    aryutt[39]=aryutt[20]*aryutt[39];
    aryutt[40]= - aryutt[68] - 5*aryutt[53] - 11./2.*aryutt[6] + 2*
      aryutt[5];
    aryutt[40]=aryutt[40]*aryutt[63];
    aryutt[47]=aryutt[61]*aryutt[7];
    aryutt[49]= - aryutt[37]*aryutt[57];
    aryutt[47]=2./3.*aryutt[47] + aryutt[49];
    aryutt[47]=MMZ*aryutt[9]*aryutt[47];
    aryutt[49]= - aryutt[62]*aryutt[51];
    aryutt[47]=aryutt[49] + aryutt[47];
    aryutt[47]=aryutt[18]*aryutt[47];
    aryutt[42]= - aryutt[52]*aryutt[42];
    aryutt[44]=aryutt[9]*aryutt[44]*aryutt[48];
    aryutt[42]=aryutt[44] + aryutt[42];
    aryutt[42]=MMZ*aryutt[42];
    aryutt[41]=MMt*aryutt[41]*aryutt[54];
    aryutt[41]=aryutt[42] - 1./6.*aryutt[41];
    aryutt[41]=aryutt[41]*pow(Pi,2);

    yuttret = aryutt[39] + aryutt[40] + 5*aryutt[41] + aryutt[43] + 
      aryutt[45] + aryutt[46] + 2*aryutt[47];
    return yuttret.real();
  }
} // namespace mr
