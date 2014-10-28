#include <bb.hpp>
std::complex<long double>
bb::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbar[29], mbbbarret;

    armbbbar[1]=double(boson);
    armbbbar[2]=pow(CW,-1);
    armbbbar[3]=pow(MMH,-1);
    armbbbar[4]=pow(MMZ,-1);
    armbbbar[5]=pow(SW,-1);
    armbbbar[6]=Tsil::I2(0,MMW,MMt,mu2);
    armbbbar[7]=Tsil::A(MMH,mu2);
    armbbbar[8]=Tsil::A(MMb,mu2);
    armbbbar[9]=pow(MMb,-1);
    armbbbar[10]=Tsil::A(MMZ,mu2);
    armbbbar[11]=Tsil::A(MMW,mu2);
    armbbbar[12]=Tsil::A(MMt,mu2);
    armbbbar[13]=Tsil::Aeps(MMW,mu2);
    armbbbar[14]=Tsil::Aeps(MMt,mu2);
    armbbbar[15]=Tsil::Aeps(MMb,mu2);
    armbbbar[16]=prot0bb0b->Tvxu(0);
    armbbbar[17]=1/(MMt - MMW);
   armbbbar[18]=pow(CW,2);
   armbbbar[19]=1 + armbbbar[18];
   armbbbar[20]=armbbbar[18]*armbbbar[19];
   armbbbar[20]=1 + armbbbar[20];
   armbbbar[20]=MMZ*armbbbar[20];
   armbbbar[18]= - 1 - armbbbar[18];
   armbbbar[21]=armbbbar[14]*armbbbar[18];
   armbbbar[22]=armbbbar[13]*armbbbar[18];
   armbbbar[23]=armbbbar[6]*armbbbar[19];
   armbbbar[24]=armbbbar[11]*armbbbar[18];
   armbbbar[18]=armbbbar[12]*armbbbar[18];
   armbbbar[18]=2*armbbbar[20] + armbbbar[18] + armbbbar[24] + 
   armbbbar[23] + armbbbar[21] + armbbbar[22];
   armbbbar[18]=MMZ*armbbbar[18];
   armbbbar[20]=armbbbar[12]*armbbbar[11];
   armbbbar[18]=armbbbar[20] + armbbbar[18];
   armbbbar[18]=MMZ*armbbbar[18];
   armbbbar[20]= - armbbbar[12]*armbbbar[11];
   armbbbar[21]= - armbbbar[6] + armbbbar[14] + armbbbar[13];
   armbbbar[22]= - 2*MMZ + armbbbar[12] + armbbbar[21] + armbbbar[11];
   armbbbar[22]=MMZ*armbbbar[22];
   armbbbar[20]=armbbbar[20] + armbbbar[22];
   armbbbar[22]=pow(armbbbar[5],2);
   armbbbar[20]=armbbbar[22]*MMZ*armbbbar[20];
   armbbbar[18]=armbbbar[18] + armbbbar[20];
   armbbbar[18]=armbbbar[17]*armbbbar[18];
   armbbbar[20]=armbbbar[6] - armbbbar[14] - armbbbar[13];
   armbbbar[23]= - armbbbar[11] + armbbbar[12];
   armbbbar[23]=armbbbar[9]*armbbbar[8]*armbbbar[23];
   armbbbar[19]=MMZ*armbbbar[19];
   armbbbar[19]=17*armbbbar[19] + 3./2.*armbbbar[23] - 11*armbbbar[12]
    + 7*armbbbar[20] - 3*armbbbar[11];
   armbbbar[19]=MMZ*armbbbar[19];
   armbbbar[20]=3*armbbbar[11];
   armbbbar[23]=armbbbar[11] - armbbbar[12];
   armbbbar[23]=armbbbar[9]*armbbbar[8]*armbbbar[23];
   armbbbar[23]= - 17*MMZ + 3./2.*armbbbar[23] + 11*armbbbar[12] + 7*
   armbbbar[21] + armbbbar[20];
   armbbbar[23]=MMZ*armbbbar[23];
   armbbbar[24]= - 4*armbbbar[11];
   armbbbar[25]=armbbbar[24] - 3./2.*armbbbar[12];
   armbbbar[25]=armbbbar[12]*armbbbar[25];
   armbbbar[23]=armbbbar[25] + armbbbar[23];
   armbbbar[23]=armbbbar[22]*armbbbar[23];
   armbbbar[18]=3*armbbbar[18] + armbbbar[19] + armbbbar[23];
   armbbbar[18]=armbbbar[17]*armbbbar[18];
   armbbbar[19]=armbbbar[9]*armbbbar[8]*armbbbar[11];
   armbbbar[23]=armbbbar[24] - 1./2.*armbbbar[12];
   armbbbar[23]=armbbbar[4]*armbbbar[12]*armbbbar[23];
   armbbbar[24]=armbbbar[9]*armbbbar[8];
   armbbbar[24]= - 10 + 1./2.*armbbbar[24];
   armbbbar[24]=MMZ*armbbbar[24];
   armbbbar[19]=3*armbbbar[24] + armbbbar[23] + 3./2.*armbbbar[19] + 9*
   armbbbar[12] + 8*armbbbar[21] + armbbbar[20];
   armbbbar[19]=armbbbar[22]*armbbbar[19];
   armbbbar[20]=pow(armbbbar[2],2);
   armbbbar[21]= - armbbbar[20]*armbbbar[11];
   armbbbar[23]= - armbbbar[12]*armbbbar[20];
   armbbbar[21]=4*armbbbar[21] + 1./2.*armbbbar[23];
   armbbbar[21]=armbbbar[4]*armbbbar[12]*armbbbar[21];
   armbbbar[23]= - armbbbar[9]*armbbbar[8];
   armbbbar[23]=10 + 1./2.*armbbbar[23];
   armbbbar[23]=MMZ*armbbbar[23];
   armbbbar[18]=armbbbar[18] + armbbbar[19] + armbbbar[21] + 3*
   armbbbar[23];
   armbbbar[18]=armbbbar[17]*armbbbar[18];
   armbbbar[19]=1./4.*MMt + 3./2.*armbbbar[7] - 2./3.*armbbbar[10];
   armbbbar[19]=armbbbar[20]*armbbbar[19];
   armbbbar[21]= - armbbbar[3]*MMt;
   armbbbar[21]=1./2. + 4*armbbbar[21];
   armbbbar[23]=armbbbar[12]*armbbbar[20]*armbbbar[21];
   armbbbar[19]=3*armbbbar[23] - 4./3.*armbbbar[10] + armbbbar[19];
   armbbbar[19]=armbbbar[9]*armbbbar[8]*armbbbar[19];
   armbbbar[23]=armbbbar[3]*pow(MMt,2);
   armbbbar[23]=32*armbbbar[23] - 103./12.*MMt - 4*armbbbar[6] - 3./2.*
   armbbbar[10] + 4*armbbbar[13] - 1./2.*armbbbar[7] + 4*armbbbar[14];
   armbbbar[24]=armbbbar[20]*armbbbar[23];
   armbbbar[25]=armbbbar[3]*MMt;
   armbbbar[25]=13./2. + 4*armbbbar[25];
   armbbbar[26]=armbbbar[20]*armbbbar[25];
   armbbbar[27]= - armbbbar[20]*armbbbar[3];
   armbbbar[28]=armbbbar[12]*armbbbar[27];
   armbbbar[26]=armbbbar[26] + 24*armbbbar[28];
   armbbbar[26]=armbbbar[12]*armbbbar[26];
   armbbbar[19]=armbbbar[19] + armbbbar[24] + armbbbar[26];
   armbbbar[19]=armbbbar[4]*armbbbar[19];
   armbbbar[24]= - armbbbar[12]*armbbbar[3];
   armbbbar[24]=armbbbar[25] + 24*armbbbar[24];
   armbbbar[24]=armbbbar[12]*armbbbar[24];
   armbbbar[25]=3*armbbbar[7] + 1./2.*MMt;
   armbbbar[21]=armbbbar[12]*armbbbar[21];
   armbbbar[21]=1./2.*armbbbar[25] + 3*armbbbar[21];
   armbbbar[21]=armbbbar[9]*armbbbar[8]*armbbbar[21];
   armbbbar[21]=armbbbar[21] + armbbbar[23] + armbbbar[24];
   armbbbar[21]=armbbbar[4]*armbbbar[21];
   armbbbar[23]= - armbbbar[10] - 2*armbbbar[11];
   armbbbar[23]=armbbbar[3]*armbbbar[23];
   armbbbar[24]=armbbbar[10] + 2*armbbbar[11];
   armbbbar[24]=armbbbar[3]*armbbbar[24];
   armbbbar[24]= - 1./4. + armbbbar[24];
   armbbbar[24]=armbbbar[9]*armbbbar[8]*armbbbar[24];
   armbbbar[25]=armbbbar[9]*armbbbar[8]*armbbbar[3];
   armbbbar[25]= - armbbbar[3] + 3*armbbbar[25];
   armbbbar[25]=MMZ*armbbbar[25];
   armbbbar[21]=2*armbbbar[25] + armbbbar[21] + 3*armbbbar[24] - 167./8.
    + armbbbar[23];
   armbbbar[21]=armbbbar[22]*armbbbar[21];
   armbbbar[22]= - armbbbar[3]*armbbbar[10];
   armbbbar[22]= - 107./216. + armbbbar[22];
   armbbbar[22]=armbbbar[20]*armbbbar[22];
   armbbbar[23]=armbbbar[3]*armbbbar[10];
   armbbbar[23]= - 31./36. + 3*armbbbar[23];
   armbbbar[23]=armbbbar[20]*armbbbar[23];
   armbbbar[23]=34./9. + armbbbar[23];
   armbbbar[23]=armbbbar[8]*armbbbar[23];
   armbbbar[24]=armbbbar[9]*pow(armbbbar[8],2);
   armbbbar[23]=32./9.*armbbbar[24] - 40./9.*armbbbar[15] + 
   armbbbar[23];
   armbbbar[23]=armbbbar[9]*armbbbar[23];
   armbbbar[24]=2*armbbbar[3] + armbbbar[27];
   armbbbar[20]=armbbbar[20]*armbbbar[3];
   armbbbar[20]= - 2*armbbbar[3] + armbbbar[20];
   armbbbar[20]=armbbbar[9]*armbbbar[8]*armbbbar[20];
   armbbbar[20]=1./3.*armbbbar[24] + armbbbar[20];
   armbbbar[20]=MMZ*armbbbar[20];
   armbbbar[24]= - 77./3. - 20*armbbbar[16];
   armbbbar[18]=armbbbar[18] + armbbbar[21] + 2*armbbbar[20] + 
   armbbbar[19] + armbbbar[23] + 2./9.*armbbbar[24] + armbbbar[22];

      mbbbarret = armbbbar[18]*armbbbar[1];
      return mbbbarret;
}
