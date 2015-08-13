#include <ZZ.hpp>
long double ZZ<OS>::x10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZbar[35], mZZbarret;

    armZZbar[1]=double(nH);
    armZZbar[2]=double(boson);
    armZZbar[3]=pow(CW,-1);
    armZZbar[4]=pow(MMZ,-1);
    armZZbar[5]=pow(SW,-1);
    armZZbar[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    armZZbar[7]=Tsil::B(MMb,MMb,MMZ,mu2);
    armZZbar[8]=Tsil::B(0,0,MMZ,mu2);
    armZZbar[9]=Tsil::A(MMt,mu2);
    armZZbar[10]=pow(MMH,-1);
    armZZbar[11]=Tsil::A(MMb,mu2);
    armZZbar[12]=double(nL + nH);
    armZZbar[13]=Tsil::B(MMZ,MMH,MMZ,mu2);
    armZZbar[14]=Tsil::B(MMW,MMW,MMZ,mu2);
    armZZbar[15]=Tsil::A(MMH,mu2);
    armZZbar[16]=Tsil::A(MMZ,mu2);
    armZZbar[17]=Tsil::A(MMW,mu2);
   armZZbar[18]=pow(armZZbar[3],2);
   armZZbar[19]=pow(armZZbar[5],2);
   armZZbar[20]=armZZbar[19] - 32./9. + 17./9.*armZZbar[18];
   armZZbar[20]=armZZbar[20]*armZZbar[1];
   armZZbar[21]=1./2.*armZZbar[19];
   armZZbar[22]=7./18.*armZZbar[18] - 32./9. - armZZbar[21];
   armZZbar[23]=armZZbar[6]*armZZbar[1];
   armZZbar[22]=armZZbar[22]*armZZbar[23];
   armZZbar[24]=armZZbar[18] + armZZbar[19];
   armZZbar[25]=armZZbar[10]*armZZbar[1];
   armZZbar[25]=6*armZZbar[25];
   armZZbar[25]=armZZbar[25]*armZZbar[24];
   armZZbar[26]= - armZZbar[9]*armZZbar[25];
   armZZbar[22]=armZZbar[26] + armZZbar[20] + armZZbar[22];
   armZZbar[22]=MMt*armZZbar[22];
   armZZbar[26]= - 1./2. - armZZbar[13];
   armZZbar[27]=MMH*armZZbar[13];
   armZZbar[27]=armZZbar[27] - armZZbar[16] + armZZbar[15];
   armZZbar[27]=armZZbar[4]*armZZbar[27];
   armZZbar[26]=1./12.*armZZbar[27] + 1./3.*armZZbar[26];
   armZZbar[26]=MMH*armZZbar[24]*armZZbar[26];
   armZZbar[27]=armZZbar[19] + 5./9.*armZZbar[18];
   armZZbar[28]= - 8./9. + armZZbar[27];
   armZZbar[28]=armZZbar[1]*armZZbar[28];
   armZZbar[25]= - MMb*armZZbar[25];
   armZZbar[25]=armZZbar[28] + armZZbar[25];
   armZZbar[25]=armZZbar[11]*armZZbar[25];
   armZZbar[28]=armZZbar[24]*armZZbar[16];
   armZZbar[29]=armZZbar[15]*armZZbar[24];
   armZZbar[27]=MMb*armZZbar[27];
   armZZbar[30]=8./9.*MMb;
   armZZbar[27]= - armZZbar[30] + armZZbar[27];
   armZZbar[27]=armZZbar[1]*armZZbar[27];
   armZZbar[20]=armZZbar[9]*armZZbar[20];
   armZZbar[31]=armZZbar[21] + 17./18.*armZZbar[18];
   armZZbar[32]= - MMb*armZZbar[31];
   armZZbar[30]= - armZZbar[30] + armZZbar[32];
   armZZbar[32]=armZZbar[7]*armZZbar[1];
   armZZbar[30]=armZZbar[30]*armZZbar[32];
   armZZbar[33]=1./6.*armZZbar[18];
   armZZbar[34]=armZZbar[33] + 4 - 5./2.*armZZbar[19];
   armZZbar[34]=armZZbar[17]*armZZbar[34];
   armZZbar[20]=armZZbar[25] + armZZbar[22] + armZZbar[34] + 
   armZZbar[30] + 1./6.*armZZbar[28] + armZZbar[20] + 1./2.*
   armZZbar[29] + armZZbar[27] + armZZbar[26];
   armZZbar[20]=armZZbar[4]*armZZbar[20];
   armZZbar[22]=pow(CW,2);
   armZZbar[22]=4*armZZbar[22];
   armZZbar[25]=armZZbar[22] + 1./12.*armZZbar[18] + 29./3. - 33./4.*
   armZZbar[19];
   armZZbar[25]=armZZbar[14]*armZZbar[25];
   armZZbar[26]=3*armZZbar[19];
   armZZbar[27]=armZZbar[18] - 2 + armZZbar[26];
   armZZbar[27]=MMZ*armZZbar[27];
   armZZbar[26]=armZZbar[17]*armZZbar[26];
   armZZbar[26]=armZZbar[26] + armZZbar[27] + 3./2.*armZZbar[28];
   armZZbar[26]=armZZbar[10]*armZZbar[26];
   armZZbar[27]= - armZZbar[33] + 8 - 59./6.*armZZbar[19];
   armZZbar[24]=armZZbar[13]*armZZbar[24];
   armZZbar[28]= - 11./9.*armZZbar[18] + 20./9. - armZZbar[19];
   armZZbar[28]=armZZbar[1]*armZZbar[8]*armZZbar[28];
   armZZbar[19]=5./3.*armZZbar[18] + armZZbar[19] - 8./3.;
   armZZbar[29]=armZZbar[8] - 1./3.;
   armZZbar[19]=armZZbar[12]*armZZbar[29]*armZZbar[19];
   armZZbar[29]= - 16./9. + armZZbar[31];
   armZZbar[23]=armZZbar[29]*armZZbar[23];
   armZZbar[18]=5./18.*armZZbar[18] - 4./9. + armZZbar[21];
   armZZbar[18]=armZZbar[18]*armZZbar[32];
   armZZbar[18]=armZZbar[20] + armZZbar[26] + armZZbar[25] + 
   armZZbar[18] + armZZbar[23] + 4./3.*armZZbar[19] + armZZbar[28] + 
   armZZbar[24] + 1./3.*armZZbar[27] + armZZbar[22];

      mZZbarret = armZZbar[18]*armZZbar[2];
      return mZZbarret.real();
}
