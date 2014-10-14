#include <ZZ.hpp>
std::complex<long double> ZZ::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myZ[24], myZret;

    myZ[1]=pow(CW,-1);
    myZ[2]=pow(MMZ,-1);
    myZ[3]=pow(SW,-1);
    myZ[4]=double(nL + nH);
    myZ[5]=std::real(Tsil::B(0,0,MMZ,mu2));
    myZ[6]=Tsil::B(MMZ,MMH,MMZ,mu2);
    myZ[7]=Tsil::B(MMW,MMW,MMZ,mu2);
    myZ[8]=Tsil::B(MMt,MMt,MMZ,mu2);
    myZ[9]=Tsil::A(MMH,mu2);
    myZ[10]=Tsil::A(MMZ,mu2);
    myZ[11]=Tsil::A(MMW,mu2);
    myZ[12]=Tsil::A(MMt,mu2);
    myZ[13]=1/( - MMW + MMH);
   myZ[14]=pow(myZ[3],2);
   myZ[15]=1./2.*myZ[14];
   myZ[16]=pow(myZ[1],2);
   myZ[17]= - myZ[15] - 32./9. + 7./18.*myZ[16];
   myZ[18]=myZ[8]*myZ[17];
   myZ[18]=myZ[18] + 41./36.*myZ[16] - 32./9. + 1./4.*myZ[14];
   myZ[18]=MMt*myZ[18];
   myZ[19]= - 1./8. - myZ[6];
   myZ[20]=MMH*myZ[6];
   myZ[20]=myZ[20] + myZ[9] - myZ[10];
   myZ[20]=myZ[2]*myZ[20];
   myZ[19]=1./12.*myZ[20] + 1./3.*myZ[19];
   myZ[20]=myZ[16] + myZ[14];
   myZ[19]=MMH*myZ[20]*myZ[19];
   myZ[21]=11./3. - 3*myZ[14];
   myZ[21]=myZ[21]*myZ[14];
   myZ[21]=myZ[21] + 11./3.*myZ[16];
   myZ[21]=myZ[10]*myZ[21];
   myZ[22]= - myZ[9]*myZ[20];
   myZ[21]=myZ[22] + myZ[21];
   myZ[17]=myZ[12]*myZ[17];
   myZ[22]= - 1 + 3./4.*myZ[14];
   myZ[22]=myZ[22]*myZ[14];
   myZ[23]=5./3.*myZ[16];
   myZ[22]=myZ[23] + 4 + myZ[22];
   myZ[22]=myZ[11]*myZ[22];
   myZ[17]=myZ[22] + myZ[17] + 1./4.*myZ[21] + myZ[19] + myZ[18];
   myZ[17]=myZ[2]*myZ[17];
   myZ[18]= - myZ[11] + myZ[9];
   myZ[19]=3./4.*myZ[13];
   myZ[18]=myZ[19]*myZ[14]*myZ[18];
   myZ[19]=pow(CW,2);
   myZ[19]=4*myZ[19];
   myZ[21]=1./12.*myZ[16] + myZ[19] + 29./3. - 33./4.*myZ[14];
   myZ[21]=myZ[7]*myZ[21];
   myZ[15]=myZ[15] - 16./9. + 17./18.*myZ[16];
   myZ[22]=myZ[23] + myZ[14] - 8./3.;
   myZ[22]=myZ[22]*myZ[4];
   myZ[23]=4./3.*myZ[22] - myZ[15];
   myZ[23]=myZ[5]*myZ[23];
   myZ[15]=myZ[8]*myZ[15];
   myZ[20]=myZ[6]*myZ[20];
   myZ[14]=8 - 209./24.*myZ[14];

      myZret = 1./3.*myZ[14] + myZ[15] + 5./72.*myZ[16] + myZ[17] + 
      myZ[18] + myZ[19] + myZ[20] + myZ[21] - 4./9.*myZ[22] + myZ[23];
      return myZret;
}
