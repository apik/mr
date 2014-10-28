#include <tt.hpp>
std::complex<long double>
tt::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armttbar[22], mttbarret;

    armttbar[1]=double(nH);
    armttbar[2]=Tsil::A(MMt,mu2);
    armttbar[3]=pow(CW,-1);
    armttbar[4]=pow(MMH,-1);
    armttbar[5]=pow(MMZ,-1);
    armttbar[6]=pow(SW,-1);
    armttbar[7]=double(boson);
    armttbar[8]=Tsil::B(MMH,MMt,MMt,mu2);
    armttbar[9]=Tsil::B(MMZ,MMt,MMt,mu2);
    armttbar[10]=pow(MMt,-1);
    armttbar[11]=Tsil::A(MMH,mu2);
    armttbar[12]=Tsil::A(MMZ,mu2);
    armttbar[13]=Tsil::A(MMW,mu2);
    armttbar[14]=std::real(Tsil::B(0,MMW,MMt,mu2));
   armttbar[15]= - armttbar[8]*MMH;
   armttbar[15]=1./2.*armttbar[15] + armttbar[11] + 1./2.*armttbar[13];
   armttbar[16]=1./4.*armttbar[14];
   armttbar[17]=armttbar[8] + armttbar[16];
   armttbar[17]=MMt*armttbar[17];
   armttbar[18]=1./2.*armttbar[2];
   armttbar[15]=armttbar[18] + 1./2.*armttbar[15] + armttbar[17];
   armttbar[15]=armttbar[5]*armttbar[15];
   armttbar[17]= - 1 - 7*armttbar[9];
   armttbar[19]=3./2.*armttbar[12] + MMZ;
   armttbar[19]=armttbar[4]*armttbar[19];
   armttbar[20]= - MMZ*armttbar[9];
   armttbar[20]=armttbar[20] - armttbar[12] + armttbar[2];
   armttbar[20]=armttbar[10]*armttbar[20];
   armttbar[17]=armttbar[15] + 17./36.*armttbar[20] + 1./36.*
   armttbar[17] + armttbar[19];
   armttbar[19]=pow(armttbar[3],2);
   armttbar[17]=armttbar[19]*armttbar[17];
   armttbar[20]= - armttbar[14] - 1./2.*armttbar[9];
   armttbar[20]=MMZ*armttbar[20];
   armttbar[18]=armttbar[20] + armttbar[18] - armttbar[13] - 1./2.*
   armttbar[12];
   armttbar[18]=armttbar[10]*armttbar[18];
   armttbar[20]=armttbar[9] - 3 + armttbar[14];
   armttbar[21]=MMZ + armttbar[13] + 1./2.*armttbar[12];
   armttbar[21]=armttbar[4]*armttbar[21];
   armttbar[15]=armttbar[15] + 1./2.*armttbar[18] + 1./4.*armttbar[20]
    + 3*armttbar[21];
   armttbar[18]=pow(armttbar[6],2);
   armttbar[15]=armttbar[18]*armttbar[15];
   armttbar[20]=armttbar[12] + 2*armttbar[2];
   armttbar[16]=armttbar[16] + 4./9.*armttbar[9];
   armttbar[16]=MMZ*armttbar[16];
   armttbar[16]=4./9.*armttbar[20] + armttbar[16];
   armttbar[16]=armttbar[10]*armttbar[16];
   armttbar[20]= - 1 + armttbar[9];
   armttbar[21]= - armttbar[4]*MMZ;
   armttbar[15]=1./2.*armttbar[15] + 1./2.*armttbar[17] + armttbar[16]
    + 8./9.*armttbar[20] + armttbar[21];
   armttbar[15]=armttbar[7]*armttbar[15];
   armttbar[16]= - armttbar[5]*armttbar[4]*armttbar[2]*MMt*armttbar[1];
   armttbar[17]=armttbar[19]*armttbar[16];
   armttbar[16]=armttbar[18]*armttbar[16];
   armttbar[16]=armttbar[17] + armttbar[16];

      mttbarret = armttbar[15] + 3*armttbar[16];
      return mttbarret;
}
