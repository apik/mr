#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZbar[27], mZZbarret;

    armZZbar[1]=double(nL + nH);
    armZZbar[2]=pow(CW,-1);
    armZZbar[3]=pow(SW,-1);
    armZZbar[4]=std::real(Tsil::B(0,0,MMZ,mu2));
    armZZbar[5]=double(nH);
    armZZbar[6]=pow(MMZ,-1);
    armZZbar[7]=Tsil::B(MMt,MMt,MMZ,mu2);
    armZZbar[8]=Tsil::A(MMt,mu2);
    armZZbar[9]=pow(MMH,-1);
    armZZbar[10]=double(nL);
    armZZbar[11]=double(boson);
    armZZbar[12]=Tsil::B(MMZ,MMH,MMZ,mu2);
    armZZbar[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    armZZbar[14]=Tsil::A(MMH,mu2);
    armZZbar[15]=Tsil::A(MMZ,mu2);
    armZZbar[16]=Tsil::A(MMW,mu2);
   armZZbar[17]=pow(armZZbar[3],2);
   armZZbar[18]=1./2.*armZZbar[17];
   armZZbar[19]=pow(armZZbar[2],2);
   armZZbar[20]=7./18.*armZZbar[19] - 32./9. - armZZbar[18];
   armZZbar[21]=armZZbar[6]*MMt;
   armZZbar[20]=armZZbar[20]*armZZbar[21];
   armZZbar[20]=armZZbar[20] + 17./18.*armZZbar[19] - 16./9. + 
   armZZbar[18];
   armZZbar[20]=armZZbar[7]*armZZbar[20];
   armZZbar[22]=11./9.*armZZbar[19] + armZZbar[17] - 20./9.;
   armZZbar[18]=5./18.*armZZbar[19] - 4./9. + armZZbar[18];
   armZZbar[18]=armZZbar[4]*armZZbar[18];
   armZZbar[23]=17./9.*armZZbar[19] + armZZbar[17] - 32./9.;
   armZZbar[24]=armZZbar[8] + MMt;
   armZZbar[23]=armZZbar[6]*armZZbar[23]*armZZbar[24];
   armZZbar[24]=armZZbar[19] + armZZbar[17];
   armZZbar[21]=armZZbar[9]*armZZbar[8]*armZZbar[24]*armZZbar[21];
   armZZbar[18]=armZZbar[20] - 6*armZZbar[21] + armZZbar[23] - 1./3.*
   armZZbar[22] + armZZbar[18];
   armZZbar[18]=armZZbar[5]*armZZbar[18];
   armZZbar[20]=1./2.*armZZbar[14];
   armZZbar[21]=armZZbar[20] - 1./2.*armZZbar[15];
   armZZbar[23]=pow(armZZbar[6],2);
   armZZbar[21]=armZZbar[23]*armZZbar[21];
   armZZbar[21]= - armZZbar[6] + armZZbar[21];
   armZZbar[21]=armZZbar[24]*armZZbar[21];
   armZZbar[25]=armZZbar[24]*armZZbar[12];
   armZZbar[23]=MMH*armZZbar[23];
   armZZbar[23]=1./4.*armZZbar[23] - armZZbar[6];
   armZZbar[23]=armZZbar[23]*armZZbar[25];
   armZZbar[21]=1./2.*armZZbar[21] + armZZbar[23];
   armZZbar[21]=MMH*armZZbar[21];
   armZZbar[23]=1./6.*armZZbar[19];
   armZZbar[21]=armZZbar[21] - armZZbar[23] + 8 - 59./6.*armZZbar[17];
   armZZbar[20]=armZZbar[20] + 1./6.*armZZbar[15];
   armZZbar[20]=armZZbar[24]*armZZbar[20];
   armZZbar[23]=armZZbar[23] + 4 - 5./2.*armZZbar[17];
   armZZbar[23]=armZZbar[16]*armZZbar[23];
   armZZbar[20]=armZZbar[23] + armZZbar[20];
   armZZbar[20]=armZZbar[6]*armZZbar[20];
   armZZbar[23]=3*armZZbar[17];
   armZZbar[26]=armZZbar[19] - 2 + armZZbar[23];
   armZZbar[26]=MMZ*armZZbar[26];
   armZZbar[24]=armZZbar[15]*armZZbar[24];
   armZZbar[23]=armZZbar[16]*armZZbar[23];
   armZZbar[23]=3./2.*armZZbar[24] + armZZbar[26] + armZZbar[23];
   armZZbar[23]=armZZbar[9]*armZZbar[23];
   armZZbar[24]=1./12.*armZZbar[19] + 29./3. - 33./4.*armZZbar[17];
   armZZbar[24]=armZZbar[13]*armZZbar[24];
   armZZbar[26]=1 + armZZbar[13];
   armZZbar[26]=armZZbar[26]*pow(CW,2);
   armZZbar[20]=armZZbar[25] + armZZbar[23] + 4*armZZbar[26] + 
   armZZbar[24] + armZZbar[20] + 1./3.*armZZbar[21];
   armZZbar[20]=armZZbar[11]*armZZbar[20];
   armZZbar[17]=armZZbar[17] - 4;
   armZZbar[17]=armZZbar[19] + 1./3.*armZZbar[17];
   armZZbar[17]=armZZbar[1]*armZZbar[17];
   armZZbar[19]=armZZbar[10]*armZZbar[22];
   armZZbar[17]=armZZbar[17] + armZZbar[19];
   armZZbar[19]=armZZbar[4] - 1./3.;
   armZZbar[17]=armZZbar[19]*armZZbar[17];

      mZZbarret = armZZbar[17] + armZZbar[18] + armZZbar[20];
      return mZZbarret;
}
