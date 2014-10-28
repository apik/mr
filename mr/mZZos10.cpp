#include <ZZ.hpp>
std::complex<long double>
ZZ<MS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZos[25], mZZosret;

    armZZos[1]=double(nL + nH);
    armZZos[2]=pow(s,-1);
    armZZos[3]=pow(c,-1);
    armZZos[4]=std::real(Tsil::B(0,0,mmZ,mu2));
    armZZos[5]=double(nH);
    armZZos[6]=pow(mmZ,-1);
    armZZos[7]=Tsil::B(mmt,mmt,mmZ,mu2);
    armZZos[8]=Tsil::A(mmt,mu2);
    armZZos[9]=pow(mmH,-1);
    armZZos[10]=double(nL);
    armZZos[11]=double(boson);
    armZZos[12]=Tsil::B(mmW,mmW,mmZ,mu2);
    armZZos[13]=Tsil::B(mmZ,mmH,mmZ,mu2);
    armZZos[14]=Tsil::A(mmW,mu2);
    armZZos[15]=Tsil::A(mmZ,mu2);
    armZZos[16]=Tsil::A(mmH,mu2);
   armZZos[17]= - mmH*armZZos[13];
   armZZos[17]=armZZos[17] - armZZos[16] + armZZos[15];
   armZZos[17]=mmH*armZZos[17];
   armZZos[18]=pow(armZZos[3],2);
   armZZos[19]=armZZos[18]*armZZos[17];
   armZZos[20]=pow(armZZos[2],2);
   armZZos[17]=armZZos[20]*armZZos[17];
   armZZos[17]=armZZos[19] + armZZos[17];
   armZZos[17]=armZZos[6]*armZZos[17];
   armZZos[19]= - 1./3.*armZZos[15];
   armZZos[21]=armZZos[19] - armZZos[16] - 1./3.*armZZos[14];
   armZZos[22]=1./2. + armZZos[13];
   armZZos[22]=1./3.*mmH*armZZos[22];
   armZZos[21]=1./2.*armZZos[21] + armZZos[22];
   armZZos[21]=armZZos[18]*armZZos[21];
   armZZos[19]=armZZos[19] - armZZos[16] + 5*armZZos[14];
   armZZos[19]=1./2.*armZZos[19] + armZZos[22];
   armZZos[19]=armZZos[20]*armZZos[19];
   armZZos[17]=1./12.*armZZos[17] + armZZos[19] - 4*armZZos[14] + 
   armZZos[21];
   armZZos[17]=armZZos[6]*armZZos[17];
   armZZos[19]=1./3. - 1./2.*armZZos[12];
   armZZos[21]= - mmZ - 3./2.*armZZos[15];
   armZZos[21]=armZZos[9]*armZZos[21];
   armZZos[19]=armZZos[21] + 1./6.*armZZos[19] - armZZos[13];
   armZZos[19]=armZZos[18]*armZZos[19];
   armZZos[21]=59./9. + 33./2.*armZZos[12];
   armZZos[22]= - 1./2.*armZZos[15] - mmZ - armZZos[14];
   armZZos[22]=armZZos[9]*armZZos[22];
   armZZos[21]=3*armZZos[22] + 1./2.*armZZos[21] - armZZos[13];
   armZZos[21]=armZZos[20]*armZZos[21];
   armZZos[22]=pow(c,2);
   armZZos[23]= - 2./3. - armZZos[22];
   armZZos[22]= - 29./3. - 4*armZZos[22];
   armZZos[22]=armZZos[12]*armZZos[22];
   armZZos[24]=armZZos[9]*mmZ;
   armZZos[17]=armZZos[17] + armZZos[21] + armZZos[19] + 2*armZZos[24]
    + 4*armZZos[23] + armZZos[22];
   armZZos[17]=armZZos[11]*armZZos[17];
   armZZos[19]= - 5./3.*armZZos[10] - armZZos[1];
   armZZos[21]=5./3.*armZZos[10] + armZZos[1];
   armZZos[21]=armZZos[4]*armZZos[21];
   armZZos[22]=armZZos[4] - 5./3. + 4*armZZos[7];
   armZZos[22]=armZZos[5]*armZZos[22];
   armZZos[19]=1./3.*armZZos[22] + 1./3.*armZZos[19] + armZZos[21];
   armZZos[21]=11./9.*armZZos[10] + armZZos[1];
   armZZos[22]= - 11./9.*armZZos[10] - armZZos[1];
   armZZos[22]=armZZos[4]*armZZos[22];
   armZZos[23]= - 5./2.*armZZos[4] + 11./3. - 17./2.*armZZos[7];
   armZZos[23]=armZZos[5]*armZZos[23];
   armZZos[21]=1./9.*armZZos[23] + 1./3.*armZZos[21] + armZZos[22];
   armZZos[21]=armZZos[18]*armZZos[21];
   armZZos[22]=armZZos[10] + 1./3.*armZZos[1];
   armZZos[23]= - armZZos[10] - 1./3.*armZZos[1];
   armZZos[23]=armZZos[4]*armZZos[23];
   armZZos[24]= - 1./2.*armZZos[4] + 1./3. - 1./2.*armZZos[7];
   armZZos[24]=armZZos[5]*armZZos[24];
   armZZos[22]=armZZos[24] + 1./3.*armZZos[22] + armZZos[23];
   armZZos[22]=armZZos[20]*armZZos[22];
   armZZos[23]= - 17 - 7./2.*armZZos[7];
   armZZos[23]=mmt*armZZos[23];
   armZZos[23]= - 17*armZZos[8] + armZZos[23];
   armZZos[24]=6*armZZos[9]*mmt*armZZos[8];
   armZZos[23]=1./9.*armZZos[23] + armZZos[24];
   armZZos[18]=armZZos[18]*armZZos[5]*armZZos[23];
   armZZos[23]= - 1 + 1./2.*armZZos[7];
   armZZos[23]=mmt*armZZos[23];
   armZZos[23]=armZZos[24] - armZZos[8] + armZZos[23];
   armZZos[20]=armZZos[20]*armZZos[5]*armZZos[23];
   armZZos[23]=1 + armZZos[7];
   armZZos[23]=mmt*armZZos[23];
   armZZos[23]=armZZos[8] + armZZos[23];
   armZZos[23]=armZZos[5]*armZZos[23];
   armZZos[18]=armZZos[20] + 32./9.*armZZos[23] + armZZos[18];
   armZZos[18]=armZZos[6]*armZZos[18];

      mZZosret = armZZos[17] + armZZos[18] + 4./3.*armZZos[19] + 
      armZZos[21] + armZZos[22];
      return mZZosret;
}
