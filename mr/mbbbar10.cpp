#include <bb.hpp>
long double bb::x10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armbbbar[27], mbbbarret;

    armbbbar[1]=double(nH);
    armbbbar[2]=double(boson);
    armbbbar[3]=Tsil::A(MMt,mu2);
    armbbbar[4]=pow(CW,-1);
    armbbbar[5]=pow(MMH,-1);
    armbbbar[6]=pow(MMZ,-1);
    armbbbar[7]=pow(SW,-1);
    armbbbar[8]=Tsil::A(MMb,mu2);
    armbbbar[9]=Tsil::B(MMH,MMb,MMb,mu2);
    armbbbar[10]=Tsil::B(MMZ,MMb,MMb,mu2);
    armbbbar[11]=pow(MMb,-1);
    armbbbar[12]=Tsil::B(MMW,MMt,MMb,mu2);
    armbbbar[13]=Tsil::A(MMH,mu2);
    armbbbar[14]=Tsil::A(MMZ,mu2);
    armbbbar[15]=Tsil::A(MMW,mu2);
   armbbbar[16]=3*armbbbar[5];
   armbbbar[17]=MMZ + armbbbar[15];
   armbbbar[17]=armbbbar[17]*armbbbar[16];
   armbbbar[18]=armbbbar[10] - 3 + armbbbar[12];
   armbbbar[19]=3./2.*armbbbar[14];
   armbbbar[20]=armbbbar[5]*armbbbar[19];
   armbbbar[17]=armbbbar[20] + armbbbar[17] + 1./4.*armbbbar[18];
   armbbbar[17]=armbbbar[2]*armbbbar[17];
   armbbbar[18]=armbbbar[10]*MMZ;
   armbbbar[18]=armbbbar[18] + armbbbar[14];
   armbbbar[20]=MMt*armbbbar[12];
   armbbbar[21]= - armbbbar[18] + armbbbar[20] + armbbbar[8];
   armbbbar[22]=armbbbar[3] - armbbbar[15];
   armbbbar[23]=MMZ*armbbbar[12];
   armbbbar[21]= - armbbbar[23] + armbbbar[22] + 1./2.*armbbbar[21];
   armbbbar[24]=armbbbar[11]*armbbbar[2];
   armbbbar[21]=armbbbar[21]*armbbbar[24];
   armbbbar[17]=armbbbar[17] + 1./2.*armbbbar[21];
   armbbbar[21]=pow(armbbbar[7],2);
   armbbbar[17]=armbbbar[17]*armbbbar[21];
   armbbbar[25]=armbbbar[8] - armbbbar[18];
   armbbbar[25]=armbbbar[25]*armbbbar[24];
   armbbbar[19]=armbbbar[19] + MMZ;
   armbbbar[19]=armbbbar[5]*armbbbar[19];
   armbbbar[19]=17./36.*armbbbar[10] - 13./36. + armbbbar[19];
   armbbbar[19]=armbbbar[2]*armbbbar[19];
   armbbbar[19]=armbbbar[19] + 5./36.*armbbbar[25];
   armbbbar[25]=pow(armbbbar[4],2);
   armbbbar[19]=armbbbar[19]*armbbbar[25];
   armbbbar[17]=armbbbar[19] + armbbbar[17];
   armbbbar[16]=armbbbar[16]*armbbbar[1];
   armbbbar[19]=armbbbar[16]*armbbbar[8];
   armbbbar[19]= - armbbbar[19] + 1./2.*armbbbar[9] + 1./8.*
   armbbbar[12];
   armbbbar[19]=armbbbar[19]*MMb;
   armbbbar[26]= - armbbbar[13] - 1./2.*armbbbar[15] + armbbbar[20] - 
   armbbbar[8];
   armbbbar[16]=armbbbar[16]*MMt;
   armbbbar[16]=armbbbar[16] - 1./8.;
   armbbbar[16]=armbbbar[16]*armbbbar[3];
   armbbbar[16]=armbbbar[19] - armbbbar[16] - 1./4.*armbbbar[26];
   armbbbar[16]=armbbbar[2]*armbbbar[16];
   armbbbar[19]=armbbbar[20] + armbbbar[22];
   armbbbar[20]=1./8.*armbbbar[24];
   armbbbar[19]=armbbbar[20]*armbbbar[19]*MMt;
   armbbbar[20]=1./8.*armbbbar[2];
   armbbbar[20]=armbbbar[9]*armbbbar[20]*MMH;
   armbbbar[16]= - armbbbar[20] + armbbbar[19] + armbbbar[16];
   armbbbar[19]=armbbbar[25] + armbbbar[21];
   armbbbar[16]=armbbbar[6]*armbbbar[16]*armbbbar[19];
   armbbbar[18]=2./9.*armbbbar[8] + 1./4.*armbbbar[23] + 1./9.*
   armbbbar[18];
   armbbbar[18]=armbbbar[18]*armbbbar[24];
   armbbbar[19]= - armbbbar[5]*MMZ;
   armbbbar[19]=2./9.*armbbbar[10] - 2./9. + armbbbar[19];
   armbbbar[19]=armbbbar[2]*armbbbar[19];

      mbbbarret = armbbbar[16] + 1./2.*armbbbar[17] + armbbbar[18] + 
      armbbbar[19];
      return mbbbarret.real();
}
