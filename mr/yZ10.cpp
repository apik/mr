#include <ZZ.hpp>
std::complex<long double> ZZ::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myZ[25], myZret;

    myZ[1]=pow(CW,-1);
    myZ[2]=pow(MMZ,-1);
    myZ[3]=pow(SW,-1);
    myZ[4]=double(nL + nH);
    myZ[5]=std::real(Tsil::B(0,0,MMZ,mu2));
    myZ[6]=double(nH);
    myZ[7]=Tsil::B(MMt,MMt,MMZ,mu2);
    myZ[8]=Tsil::A(MMt,mu2);
    myZ[9]=double(nL);
    myZ[10]=Tsil::B(MMZ,MMH,MMZ,mu2);
    myZ[11]=Tsil::B(MMW,MMW,MMZ,mu2);
    myZ[12]=Tsil::A(MMH,mu2);
    myZ[13]=Tsil::A(MMZ,mu2);
    myZ[14]=Tsil::A(MMW,mu2);
    myZ[15]=1/( - MMW + MMH);
   myZ[16]= - 1./8. - myZ[10];
   myZ[17]=MMH*myZ[10];
   myZ[17]= - myZ[13] + myZ[12] + myZ[17];
   myZ[17]=myZ[2]*myZ[17];
   myZ[16]=1./12.*myZ[17] + 1./3.*myZ[16];
   myZ[17]=pow(myZ[1],2);
   myZ[18]=pow(myZ[3],2);
   myZ[19]=myZ[17] + myZ[18];
   myZ[16]=myZ[16]*myZ[19]*MMH;
   myZ[20]=myZ[12]*myZ[19];
   myZ[21]=11./3. - 3*myZ[18];
   myZ[21]=myZ[21]*myZ[18];
   myZ[21]=myZ[21] + 11./3.*myZ[17];
   myZ[21]=myZ[13]*myZ[21];
   myZ[20]=myZ[20] - myZ[21];
   myZ[21]=1./2.*myZ[18];
   myZ[22]= - myZ[21] - 32./9. + 7./18.*myZ[17];
   myZ[23]=myZ[7]*myZ[22];
   myZ[23]=myZ[23] + 41./36.*myZ[17] - 32./9. + 1./4.*myZ[18];
   myZ[23]=MMt*myZ[23];
   myZ[22]=myZ[8]*myZ[22];
   myZ[22]=myZ[22] + myZ[23];
   myZ[22]=myZ[6]*myZ[22];
   myZ[23]= - 1 + 3./4.*myZ[18];
   myZ[23]=myZ[23]*myZ[18];
   myZ[23]=5./3.*myZ[17] + 4 + myZ[23];
   myZ[23]=myZ[14]*myZ[23];
   myZ[16]=myZ[23] + myZ[22] - 1./4.*myZ[20] + myZ[16];
   myZ[16]=myZ[2]*myZ[16];
   myZ[20]=11./9.*myZ[17] + myZ[18] - 20./9.;
   myZ[22]=17./18.*myZ[17] - 16./9. + myZ[21];
   myZ[22]=myZ[7]*myZ[22];
   myZ[21]=5./18.*myZ[17] - 4./9. + myZ[21];
   myZ[21]=myZ[5]*myZ[21];
   myZ[21]=myZ[21] - 1./3.*myZ[20] + myZ[22];
   myZ[21]=myZ[6]*myZ[21];
   myZ[20]=myZ[9]*myZ[20];
   myZ[22]=myZ[18] - 4;
   myZ[22]=myZ[17] + 1./3.*myZ[22];
   myZ[22]=myZ[4]*myZ[22];
   myZ[20]=myZ[22] + myZ[20];
   myZ[22]=myZ[5] - 1./3.;
   myZ[20]=myZ[22]*myZ[20];
   myZ[22]= - myZ[14] + myZ[12];
   myZ[23]=3./4.*myZ[15];
   myZ[22]=myZ[23]*myZ[18]*myZ[22];
   myZ[23]=pow(CW,2);
   myZ[23]=4*myZ[23];
   myZ[24]=myZ[23] + 1./12.*myZ[17] + 29./3. - 33./4.*myZ[18];
   myZ[24]=myZ[11]*myZ[24];
   myZ[17]=5./24.*myZ[17] + 8 - 209./24.*myZ[18];
   myZ[18]=myZ[10]*myZ[19];

      myZret = myZ[16] + 1./3.*myZ[17] + myZ[18] + myZ[20] + myZ[21] + 
      myZ[22] + myZ[23] + myZ[24];
      return myZret;
}
