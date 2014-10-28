#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZbar[25], mZZbarret;

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
   armZZbar[17]=MMH*armZZbar[12];
   armZZbar[17]=armZZbar[17] + armZZbar[14] - armZZbar[15];
   armZZbar[17]=MMH*armZZbar[17];
   armZZbar[18]=pow(armZZbar[2],2);
   armZZbar[19]=armZZbar[18]*armZZbar[17];
   armZZbar[20]=pow(armZZbar[3],2);
   armZZbar[17]=armZZbar[20]*armZZbar[17];
   armZZbar[17]=armZZbar[19] + armZZbar[17];
   armZZbar[17]=armZZbar[6]*armZZbar[17];
   armZZbar[19]=1./3.*armZZbar[15];
   armZZbar[21]=armZZbar[19] + 1./3.*armZZbar[16] + armZZbar[14];
   armZZbar[22]= - 1./2. - armZZbar[12];
   armZZbar[22]=1./3.*MMH*armZZbar[22];
   armZZbar[21]=1./2.*armZZbar[21] + armZZbar[22];
   armZZbar[21]=armZZbar[18]*armZZbar[21];
   armZZbar[19]=armZZbar[19] - 5*armZZbar[16] + armZZbar[14];
   armZZbar[19]=1./2.*armZZbar[19] + armZZbar[22];
   armZZbar[19]=armZZbar[20]*armZZbar[19];
   armZZbar[17]=1./12.*armZZbar[17] + armZZbar[19] + 4*armZZbar[16] + 
   armZZbar[21];
   armZZbar[17]=armZZbar[6]*armZZbar[17];
   armZZbar[19]= - 1./3. + 1./2.*armZZbar[13];
   armZZbar[21]=MMZ + 3./2.*armZZbar[15];
   armZZbar[21]=armZZbar[9]*armZZbar[21];
   armZZbar[19]=armZZbar[21] + 1./6.*armZZbar[19] + armZZbar[12];
   armZZbar[19]=armZZbar[18]*armZZbar[19];
   armZZbar[21]= - 59./9. - 33./2.*armZZbar[13];
   armZZbar[22]=1./2.*armZZbar[15] + MMZ + armZZbar[16];
   armZZbar[22]=armZZbar[9]*armZZbar[22];
   armZZbar[21]=3*armZZbar[22] + 1./2.*armZZbar[21] + armZZbar[12];
   armZZbar[21]=armZZbar[20]*armZZbar[21];
   armZZbar[22]=pow(CW,2);
   armZZbar[23]=2./3. + armZZbar[22];
   armZZbar[22]=29./3. + 4*armZZbar[22];
   armZZbar[22]=armZZbar[13]*armZZbar[22];
   armZZbar[24]= - armZZbar[9]*MMZ;
   armZZbar[17]=armZZbar[17] + armZZbar[21] + armZZbar[19] + 2*
   armZZbar[24] + 4*armZZbar[23] + armZZbar[22];
   armZZbar[17]=armZZbar[11]*armZZbar[17];
   armZZbar[19]=5./3.*armZZbar[10] + armZZbar[1];
   armZZbar[21]= - 5./3.*armZZbar[10] - armZZbar[1];
   armZZbar[21]=armZZbar[4]*armZZbar[21];
   armZZbar[22]= - armZZbar[4] + 5./3. - 4*armZZbar[7];
   armZZbar[22]=armZZbar[5]*armZZbar[22];
   armZZbar[19]=1./3.*armZZbar[22] + 1./3.*armZZbar[19] + armZZbar[21];
   armZZbar[21]= - 11./9.*armZZbar[10] - armZZbar[1];
   armZZbar[22]=11./9.*armZZbar[10] + armZZbar[1];
   armZZbar[22]=armZZbar[4]*armZZbar[22];
   armZZbar[23]=5./2.*armZZbar[4] - 11./3. + 17./2.*armZZbar[7];
   armZZbar[23]=armZZbar[5]*armZZbar[23];
   armZZbar[21]=1./9.*armZZbar[23] + 1./3.*armZZbar[21] + armZZbar[22];
   armZZbar[21]=armZZbar[18]*armZZbar[21];
   armZZbar[22]= - armZZbar[10] - 1./3.*armZZbar[1];
   armZZbar[23]=armZZbar[10] + 1./3.*armZZbar[1];
   armZZbar[23]=armZZbar[4]*armZZbar[23];
   armZZbar[24]=1./2.*armZZbar[4] - 1./3. + 1./2.*armZZbar[7];
   armZZbar[24]=armZZbar[5]*armZZbar[24];
   armZZbar[22]=armZZbar[24] + 1./3.*armZZbar[22] + armZZbar[23];
   armZZbar[22]=armZZbar[20]*armZZbar[22];
   armZZbar[23]=17 + 7./2.*armZZbar[7];
   armZZbar[23]=MMt*armZZbar[23];
   armZZbar[23]=17*armZZbar[8] + armZZbar[23];
   armZZbar[24]= - 6*armZZbar[9]*MMt*armZZbar[8];
   armZZbar[23]=1./9.*armZZbar[23] + armZZbar[24];
   armZZbar[18]=armZZbar[18]*armZZbar[5]*armZZbar[23];
   armZZbar[23]=1 - 1./2.*armZZbar[7];
   armZZbar[23]=MMt*armZZbar[23];
   armZZbar[23]=armZZbar[24] + armZZbar[8] + armZZbar[23];
   armZZbar[20]=armZZbar[20]*armZZbar[5]*armZZbar[23];
   armZZbar[23]= - 1 - armZZbar[7];
   armZZbar[23]=MMt*armZZbar[23];
   armZZbar[23]= - armZZbar[8] + armZZbar[23];
   armZZbar[23]=armZZbar[5]*armZZbar[23];
   armZZbar[18]=armZZbar[20] + 32./9.*armZZbar[23] + armZZbar[18];
   armZZbar[18]=armZZbar[6]*armZZbar[18];

      mZZbarret = armZZbar[17] + armZZbar[18] + 4./3.*armZZbar[19] + 
      armZZbar[21] + armZZbar[22];
      return mZZbarret;
}
