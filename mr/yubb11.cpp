#include <bb.hpp>
std::complex<long double>
bb::my11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryubb[29], yubbret;

    aryubb[1]=double(boson);
    aryubb[2]=pow(CW,-1);
    aryubb[3]=pow(MMH,-1);
    aryubb[4]=pow(MMZ,-1);
    aryubb[5]=pow(SW,-1);
    aryubb[6]=Tsil::I2(0,MMW,MMt,mu2);
    aryubb[7]=Tsil::A(MMH,mu2);
    aryubb[8]=Tsil::A(MMb,mu2);
    aryubb[9]=pow(MMb,-1);
    aryubb[10]=Tsil::A(MMZ,mu2);
    aryubb[11]=Tsil::A(MMW,mu2);
    aryubb[12]=Tsil::A(MMt,mu2);
    aryubb[13]=Tsil::Aeps(MMW,mu2);
    aryubb[14]=Tsil::Aeps(MMt,mu2);
    aryubb[15]=Tsil::Aeps(MMb,mu2);
    aryubb[16]=prot0bb0b->Tvxu(0);
    aryubb[17]=1/(MMt - MMW);
   aryubb[18]=pow(CW,2);
   aryubb[19]=1 + aryubb[18];
   aryubb[20]=aryubb[18]*aryubb[19];
   aryubb[20]=1 + aryubb[20];
   aryubb[20]=MMZ*aryubb[20];
   aryubb[18]= - 1 - aryubb[18];
   aryubb[21]=aryubb[14]*aryubb[18];
   aryubb[22]=aryubb[13]*aryubb[18];
   aryubb[23]=aryubb[6]*aryubb[19];
   aryubb[24]=aryubb[11]*aryubb[18];
   aryubb[18]=aryubb[12]*aryubb[18];
   aryubb[18]=2*aryubb[20] + aryubb[18] + aryubb[24] + aryubb[23] + 
   aryubb[21] + aryubb[22];
   aryubb[18]=MMZ*aryubb[18];
   aryubb[20]=aryubb[12]*aryubb[11];
   aryubb[18]=aryubb[20] + aryubb[18];
   aryubb[18]=MMZ*aryubb[18];
   aryubb[20]= - aryubb[12]*aryubb[11];
   aryubb[21]= - aryubb[6] + aryubb[14] + aryubb[13];
   aryubb[22]= - 2*MMZ + aryubb[12] + aryubb[21] + aryubb[11];
   aryubb[22]=MMZ*aryubb[22];
   aryubb[20]=aryubb[20] + aryubb[22];
   aryubb[22]=pow(aryubb[5],2);
   aryubb[20]=aryubb[22]*MMZ*aryubb[20];
   aryubb[18]=aryubb[18] + aryubb[20];
   aryubb[18]=aryubb[17]*aryubb[18];
   aryubb[20]=aryubb[6] - aryubb[14] - aryubb[13];
   aryubb[23]= - aryubb[11] + aryubb[12];
   aryubb[23]=aryubb[9]*aryubb[8]*aryubb[23];
   aryubb[19]=MMZ*aryubb[19];
   aryubb[19]=17*aryubb[19] + 3./2.*aryubb[23] - 11*aryubb[12] + 7*
   aryubb[20] - 3*aryubb[11];
   aryubb[19]=MMZ*aryubb[19];
   aryubb[20]=3*aryubb[11];
   aryubb[23]=aryubb[11] - aryubb[12];
   aryubb[23]=aryubb[9]*aryubb[8]*aryubb[23];
   aryubb[23]= - 17*MMZ + 3./2.*aryubb[23] + 11*aryubb[12] + 7*
   aryubb[21] + aryubb[20];
   aryubb[23]=MMZ*aryubb[23];
   aryubb[24]= - 4*aryubb[11];
   aryubb[25]=aryubb[24] - 3./2.*aryubb[12];
   aryubb[25]=aryubb[12]*aryubb[25];
   aryubb[23]=aryubb[25] + aryubb[23];
   aryubb[23]=aryubb[22]*aryubb[23];
   aryubb[18]=3*aryubb[18] + aryubb[19] + aryubb[23];
   aryubb[18]=aryubb[17]*aryubb[18];
   aryubb[19]=aryubb[9]*aryubb[8]*aryubb[11];
   aryubb[23]=aryubb[24] - 1./2.*aryubb[12];
   aryubb[23]=aryubb[4]*aryubb[12]*aryubb[23];
   aryubb[24]=aryubb[9]*aryubb[8];
   aryubb[24]= - 10 + 1./2.*aryubb[24];
   aryubb[24]=MMZ*aryubb[24];
   aryubb[19]=3*aryubb[24] + aryubb[23] + 3./2.*aryubb[19] + 9*
   aryubb[12] + 8*aryubb[21] + aryubb[20];
   aryubb[19]=aryubb[22]*aryubb[19];
   aryubb[20]=pow(aryubb[2],2);
   aryubb[21]= - aryubb[20]*aryubb[11];
   aryubb[23]= - aryubb[12]*aryubb[20];
   aryubb[21]=4*aryubb[21] + 1./2.*aryubb[23];
   aryubb[21]=aryubb[4]*aryubb[12]*aryubb[21];
   aryubb[23]= - aryubb[9]*aryubb[8];
   aryubb[23]=10 + 1./2.*aryubb[23];
   aryubb[23]=MMZ*aryubb[23];
   aryubb[18]=aryubb[18] + aryubb[19] + aryubb[21] + 3*aryubb[23];
   aryubb[18]=aryubb[17]*aryubb[18];
   aryubb[19]=1./4.*MMt + 3./2.*aryubb[7] - 2./3.*aryubb[10];
   aryubb[19]=aryubb[20]*aryubb[19];
   aryubb[21]= - aryubb[3]*MMt;
   aryubb[21]=1./2. + 4*aryubb[21];
   aryubb[23]=aryubb[12]*aryubb[20]*aryubb[21];
   aryubb[19]=3*aryubb[23] - 4./3.*aryubb[10] + aryubb[19];
   aryubb[19]=aryubb[9]*aryubb[8]*aryubb[19];
   aryubb[23]=aryubb[3]*pow(MMt,2);
   aryubb[23]=32*aryubb[23] - 103./12.*MMt - 4*aryubb[6] - 3./2.*
   aryubb[10] + 4*aryubb[13] - 1./2.*aryubb[7] + 4*aryubb[14];
   aryubb[24]=aryubb[20]*aryubb[23];
   aryubb[25]=aryubb[3]*MMt;
   aryubb[25]=13./2. + 4*aryubb[25];
   aryubb[26]=aryubb[20]*aryubb[25];
   aryubb[27]= - aryubb[20]*aryubb[3];
   aryubb[28]=aryubb[12]*aryubb[27];
   aryubb[26]=aryubb[26] + 24*aryubb[28];
   aryubb[26]=aryubb[12]*aryubb[26];
   aryubb[19]=aryubb[19] + aryubb[24] + aryubb[26];
   aryubb[19]=aryubb[4]*aryubb[19];
   aryubb[24]= - aryubb[12]*aryubb[3];
   aryubb[24]=aryubb[25] + 24*aryubb[24];
   aryubb[24]=aryubb[12]*aryubb[24];
   aryubb[25]=3*aryubb[7] + 1./2.*MMt;
   aryubb[21]=aryubb[12]*aryubb[21];
   aryubb[21]=1./2.*aryubb[25] + 3*aryubb[21];
   aryubb[21]=aryubb[9]*aryubb[8]*aryubb[21];
   aryubb[21]=aryubb[21] + aryubb[23] + aryubb[24];
   aryubb[21]=aryubb[4]*aryubb[21];
   aryubb[23]= - aryubb[10] - 2*aryubb[11];
   aryubb[23]=aryubb[3]*aryubb[23];
   aryubb[24]=aryubb[10] + 2*aryubb[11];
   aryubb[24]=aryubb[3]*aryubb[24];
   aryubb[24]= - 1./4. + aryubb[24];
   aryubb[24]=aryubb[9]*aryubb[8]*aryubb[24];
   aryubb[25]=aryubb[9]*aryubb[8]*aryubb[3];
   aryubb[25]= - aryubb[3] + 3*aryubb[25];
   aryubb[25]=MMZ*aryubb[25];
   aryubb[21]=2*aryubb[25] + aryubb[21] + 3*aryubb[24] - 167./8. + 
   aryubb[23];
   aryubb[21]=aryubb[22]*aryubb[21];
   aryubb[22]= - aryubb[3]*aryubb[10];
   aryubb[22]= - 107./216. + aryubb[22];
   aryubb[22]=aryubb[20]*aryubb[22];
   aryubb[23]=aryubb[3]*aryubb[10];
   aryubb[23]= - 31./36. + 3*aryubb[23];
   aryubb[23]=aryubb[20]*aryubb[23];
   aryubb[23]=34./9. + aryubb[23];
   aryubb[23]=aryubb[8]*aryubb[23];
   aryubb[24]=aryubb[9]*pow(aryubb[8],2);
   aryubb[23]=32./9.*aryubb[24] - 40./9.*aryubb[15] + aryubb[23];
   aryubb[23]=aryubb[9]*aryubb[23];
   aryubb[24]=2*aryubb[3] + aryubb[27];
   aryubb[20]=aryubb[20]*aryubb[3];
   aryubb[20]= - 2*aryubb[3] + aryubb[20];
   aryubb[20]=aryubb[9]*aryubb[8]*aryubb[20];
   aryubb[20]=1./3.*aryubb[24] + aryubb[20];
   aryubb[20]=MMZ*aryubb[20];
   aryubb[24]= - 77./3. - 20*aryubb[16];
   aryubb[18]=aryubb[18] + aryubb[21] + 2*aryubb[20] + aryubb[19] + 
   aryubb[23] + 2./9.*aryubb[24] + aryubb[22];

      yubbret = aryubb[18]*aryubb[1];
      return yubbret;
}
