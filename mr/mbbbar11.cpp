#include <bb.hpp>
long double bb::x11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbar[37], mbbbarret;

    armbbbar[1]=double(nH);
    armbbbar[2]=pow(CW,-1);
    armbbbar[3]=pow(MMH,-1);
    armbbbar[4]=pow(MMZ,-1);
    armbbbar[5]=pow(SW,-1);
    armbbbar[6]=Tsil::A(MMt,mu2);
    armbbbar[7]=Tsil::A(MMb,mu2);
    armbbbar[8]=pow(MMb,-1);
    armbbbar[9]=double(boson);
    armbbbar[10]=Tsil::I2(0,MMW,MMt,mu2);
    armbbbar[11]=Tsil::A(MMH,mu2);
    armbbbar[12]=Tsil::A(MMZ,mu2);
    armbbbar[13]=Tsil::A(MMW,mu2);
    armbbbar[14]=Tsil::Aeps(MMW,mu2);
    armbbbar[15]=Tsil::Aeps(MMt,mu2);
    armbbbar[16]=Tsil::Aeps(MMb,mu2);
    armbbbar[17]=prot0bb0b->Suxv(0);
    armbbbar[18]=1/(MMt - MMW);
   armbbbar[19]= - armbbbar[15] + armbbbar[10] - armbbbar[14];
   armbbbar[20]=armbbbar[8]*armbbbar[7];
   armbbbar[21]=3./2.*armbbbar[20];
   armbbbar[22]=11 - armbbbar[21];
   armbbbar[22]=armbbbar[6]*armbbbar[22];
   armbbbar[22]=armbbbar[22] - 7*armbbbar[19];
   armbbbar[22]=armbbbar[22]*armbbbar[18];
   armbbbar[23]=1./2.*armbbbar[20];
   armbbbar[24]=armbbbar[23] - 10;
   armbbbar[22]=armbbbar[22] + 3*armbbbar[24];
   armbbbar[22]=armbbbar[22]*armbbbar[18];
   armbbbar[24]=3*armbbbar[20];
   armbbbar[25]=armbbbar[24] - 1;
   armbbbar[26]=2*armbbbar[3];
   armbbbar[27]=armbbbar[25]*armbbbar[26];
   armbbbar[28]=armbbbar[6] - armbbbar[19];
   armbbbar[29]=3*armbbbar[18];
   armbbbar[28]=armbbbar[28]*armbbbar[29];
   armbbbar[28]=armbbbar[28] - 17;
   armbbbar[29]=pow(armbbbar[18],2);
   armbbbar[28]=armbbbar[28]*armbbbar[29];
   armbbbar[30]=pow(armbbbar[18],3);
   armbbbar[31]=armbbbar[30]*MMZ;
   armbbbar[32]=armbbbar[28] - 6*armbbbar[31];
   armbbbar[32]=MMZ*armbbbar[32];
   armbbbar[32]=armbbbar[32] + armbbbar[27] + armbbbar[22];
   armbbbar[32]=MMZ*armbbbar[32];
   armbbbar[33]=armbbbar[20] - 103./3.;
   armbbbar[34]=1./4.*MMt;
   armbbbar[33]=armbbbar[33]*armbbbar[34];
   armbbbar[34]=armbbbar[11] + armbbbar[6];
   armbbbar[21]=armbbbar[34]*armbbbar[21];
   armbbbar[34]=pow(armbbbar[6],2);
   armbbbar[35]=armbbbar[34]*armbbbar[18];
   armbbbar[36]=armbbbar[35] + armbbbar[11];
   armbbbar[21]=armbbbar[33] + armbbbar[21] + 13./2.*armbbbar[6] - 1./2.
   *armbbbar[36] - 4*armbbbar[19];
   armbbbar[33]= - 3./2.*armbbbar[12] + armbbbar[21];
   armbbbar[33]=armbbbar[4]*armbbbar[33];
   armbbbar[19]= - 3./2.*armbbbar[35] + 9*armbbbar[6] - 8*armbbbar[19];
   armbbbar[19]=armbbbar[18]*armbbbar[19];
   armbbbar[35]=armbbbar[25]*armbbbar[12]*armbbbar[3];
   armbbbar[24]= - 167./2. - armbbbar[24];
   armbbbar[19]=armbbbar[33] + armbbbar[32] + armbbbar[19] + 1./4.*
   armbbbar[24] + armbbbar[35];
   armbbbar[24]=pow(armbbbar[5],2);
   armbbbar[19]=armbbbar[19]*armbbbar[24];
   armbbbar[32]=pow(CW,2);
   armbbbar[33]=armbbbar[32] + 1;
   armbbbar[28]= - armbbbar[33]*armbbbar[28];
   armbbbar[33]=armbbbar[30]*armbbbar[33];
   armbbbar[32]=armbbbar[33]*armbbbar[32];
   armbbbar[30]=armbbbar[30] + armbbbar[32];
   armbbbar[30]=MMZ*armbbbar[30];
   armbbbar[28]=6*armbbbar[30] + armbbbar[28];
   armbbbar[28]=MMZ*armbbbar[28];
   armbbbar[30]=pow(armbbbar[2],2);
   armbbbar[26]=armbbbar[30]*armbbbar[26];
   armbbbar[26]=armbbbar[26] - 4*armbbbar[3];
   armbbbar[32]=armbbbar[20] - 1./3.;
   armbbbar[26]=armbbbar[32]*armbbbar[26];
   armbbbar[22]=armbbbar[28] - armbbbar[22] + armbbbar[26];
   armbbbar[22]=MMZ*armbbbar[22];
   armbbbar[26]= - 3./2. - 2./3.*armbbbar[20];
   armbbbar[26]=armbbbar[12]*armbbbar[26];
   armbbbar[21]=armbbbar[26] + armbbbar[21];
   armbbbar[21]=armbbbar[21]*armbbbar[30];
   armbbbar[26]=armbbbar[12]*armbbbar[20];
   armbbbar[21]= - 4./3.*armbbbar[26] + armbbbar[21];
   armbbbar[21]=armbbbar[4]*armbbbar[21];
   armbbbar[23]=armbbbar[23] + 1;
   armbbbar[26]=armbbbar[18]*armbbbar[6];
   armbbbar[28]=4*armbbbar[26];
   armbbbar[32]=3*armbbbar[23] - armbbbar[28];
   armbbbar[32]=armbbbar[18]*armbbbar[32];
   armbbbar[23]=armbbbar[23] - armbbbar[26];
   armbbbar[23]=armbbbar[23]*armbbbar[29];
   armbbbar[26]=armbbbar[23] + armbbbar[31];
   armbbbar[29]=3*MMZ;
   armbbbar[26]=armbbbar[26]*armbbbar[29];
   armbbbar[28]=armbbbar[28]*armbbbar[4];
   armbbbar[26]= - armbbbar[28] + armbbbar[26] + armbbbar[27] + 
   armbbbar[32];
   armbbbar[26]=armbbbar[26]*armbbbar[24];
   armbbbar[27]= - MMZ*armbbbar[33];
   armbbbar[23]= - armbbbar[23] + armbbbar[27];
   armbbbar[23]=armbbbar[23]*armbbbar[29];
   armbbbar[27]= - armbbbar[30]*armbbbar[28];
   armbbbar[23]=armbbbar[26] + armbbbar[23] + armbbbar[27];
   armbbbar[23]=armbbbar[13]*armbbbar[23];
   armbbbar[20]= - 107./6. - 31*armbbbar[20];
   armbbbar[20]=1./36.*armbbbar[20] + armbbbar[35];
   armbbbar[20]=armbbbar[20]*armbbbar[30];
   armbbbar[26]=armbbbar[8]*pow(armbbbar[7],2);
   armbbbar[26]=16*armbbbar[26] - 13*armbbbar[7] - 20*armbbbar[16];
   armbbbar[26]=armbbbar[8]*armbbbar[26];
   armbbbar[26]= - 109./3. + 2*armbbbar[26];
   armbbbar[27]=armbbbar[17]*armbbbar[8];
   armbbbar[19]=armbbbar[23] + 40./9.*armbbbar[27] + armbbbar[19] + 
   armbbbar[21] + armbbbar[22] + 1./9.*armbbbar[26] + armbbbar[20];
   armbbbar[19]=armbbbar[9]*armbbbar[19];
   armbbbar[20]=armbbbar[25]*armbbbar[6];
   armbbbar[20]= - armbbbar[20] + 8*MMt;
   armbbbar[20]=MMt*armbbbar[20];
   armbbbar[20]=armbbbar[20] - 6*armbbbar[34];
   armbbbar[21]=armbbbar[30] + armbbbar[24];
   armbbbar[20]=armbbbar[1]*armbbbar[4]*armbbbar[21]*armbbbar[3]*
   armbbbar[20];

      mbbbarret = armbbbar[19] + 4*armbbbar[20];
      return mbbbarret.real();
}
