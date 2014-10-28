#include <HH.hpp>
std::complex<long double>
HH<OS>::m10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armHHbar[21], mHHbarret;

    armHHbar[1]=double(nH);
    armHHbar[2]=Tsil::B(MMt,MMt,MMH,mu2);
    armHHbar[3]=pow(CW,-1);
    armHHbar[4]=pow(MMH,-1);
    armHHbar[5]=pow(MMZ,-1);
    armHHbar[6]=pow(SW,-1);
    armHHbar[7]=Tsil::A(MMt,mu2);
    armHHbar[8]=double(boson);
    armHHbar[9]=Tsil::B(MMH,MMH,MMH,mu2);
    armHHbar[10]=Tsil::B(MMZ,MMZ,MMH,mu2);
    armHHbar[11]=Tsil::B(MMW,MMW,MMH,mu2);
    armHHbar[12]=Tsil::A(MMH,mu2);
    armHHbar[13]=Tsil::A(MMZ,mu2);
    armHHbar[14]=Tsil::A(MMW,mu2);
   armHHbar[15]=1./2.*armHHbar[13];
   armHHbar[16]=9./2.*armHHbar[9] + armHHbar[11];
   armHHbar[16]=MMH*armHHbar[16];
   armHHbar[17]=armHHbar[10]*MMH;
   armHHbar[16]=1./4.*armHHbar[17] + 1./2.*armHHbar[16] + armHHbar[15]
    + 3./2.*armHHbar[12] + armHHbar[14];
   armHHbar[17]=pow(armHHbar[3],2);
   armHHbar[18]=armHHbar[17]*armHHbar[16];
   armHHbar[19]=pow(armHHbar[6],2);
   armHHbar[16]=armHHbar[19]*armHHbar[16];
   armHHbar[16]=armHHbar[18] + armHHbar[16];
   armHHbar[16]=armHHbar[5]*armHHbar[16];
   armHHbar[18]=armHHbar[11]*MMZ;
   armHHbar[20]=armHHbar[10]*MMZ;
   armHHbar[15]=1./2.*armHHbar[20] + armHHbar[18] + armHHbar[14] + 
   armHHbar[15];
   armHHbar[15]=armHHbar[4]*armHHbar[15];
   armHHbar[15]=3*armHHbar[15] - armHHbar[11] - 1./2.*armHHbar[10];
   armHHbar[15]=armHHbar[19]*armHHbar[15];
   armHHbar[18]=armHHbar[13] + armHHbar[20];
   armHHbar[18]=armHHbar[4]*armHHbar[18];
   armHHbar[18]= - armHHbar[10] + 3*armHHbar[18];
   armHHbar[18]=armHHbar[17]*armHHbar[18];
   armHHbar[20]= - armHHbar[4]*armHHbar[11]*MMZ;
   armHHbar[15]=1./2.*armHHbar[16] + armHHbar[15] + 3*armHHbar[20] + 1./
   2.*armHHbar[18];
   armHHbar[15]=armHHbar[8]*armHHbar[15];
   armHHbar[16]=armHHbar[1]*MMt*armHHbar[2];
   armHHbar[18]= - MMt*armHHbar[2];
   armHHbar[18]= - armHHbar[7] + armHHbar[18];
   armHHbar[18]=armHHbar[4]*armHHbar[1]*MMt*armHHbar[18];
   armHHbar[16]=1./2.*armHHbar[16] + 2*armHHbar[18];
   armHHbar[17]=armHHbar[17]*armHHbar[16];
   armHHbar[16]=armHHbar[19]*armHHbar[16];
   armHHbar[16]=armHHbar[17] + armHHbar[16];
   armHHbar[16]=armHHbar[5]*armHHbar[16];

      mHHbarret = armHHbar[15] + 3*armHHbar[16];
      return mHHbarret;
}
