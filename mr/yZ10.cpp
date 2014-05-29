#include <ZZ.hpp>
std::complex<long double> ZZ::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myZ[29];

    myZ[1]=pow(CW,-1);
    myZ[2]=pow(MMZ,-1);
    myZ[3]=pow(SW,-1);
    myZ[4]=double(nL + nH);
    myZ[5]=Tsil::B(0,0,MMZ,mu2);
    myZ[6]=double(nH);
    myZ[7]=Tsil::B(MMt,MMt,MMZ,mu2);
    myZ[8]=Tsil::B(MMb,MMb,MMZ,mu2);
    myZ[9]=Tsil::A(MMt,mu2);
    myZ[10]=Tsil::A(MMb,mu2);
    myZ[11]=Tsil::B(MMZ,MMH,MMZ,mu2);
    myZ[12]=Tsil::B(MMW,MMW,MMZ,mu2);
    myZ[13]=Tsil::A(MMH,mu2);
    myZ[14]=Tsil::A(MMZ,mu2);
    myZ[15]=Tsil::A(MMW,mu2);
    myZ[16]=1/( - MMb + MMt);
    myZ[17]=1/( - MMW + MMH);
   myZ[18]=pow(myZ[1],2);
   myZ[19]=pow(myZ[3],2);
   myZ[20]=5./9.*myZ[18] + myZ[19] - 8./9.;
   myZ[21]=myZ[18] + myZ[19];
   myZ[22]=MMb*myZ[16];
   myZ[23]=3./4.*myZ[22];
   myZ[23]=myZ[23]*myZ[21];
   myZ[24]=1./2.*myZ[19];
   myZ[25]=myZ[24] + 17./18.*myZ[18];
   myZ[26]= - 8./9. - myZ[25];
   myZ[26]=myZ[8]*myZ[26];
   myZ[26]=myZ[26] + myZ[23] + myZ[20];
   myZ[26]=MMb*myZ[26];
   myZ[27]= - myZ[24] - 32./9. + 7./18.*myZ[18];
   myZ[28]=myZ[7]*myZ[27];
   myZ[23]=myZ[28] - myZ[23] + 41./36.*myZ[18] - 32./9. + 1./4.*myZ[19]
   ;
   myZ[23]=MMt*myZ[23];
   myZ[22]=3./2.*myZ[22];
   myZ[22]=myZ[22]*myZ[21];
   myZ[27]= - myZ[22] + myZ[27];
   myZ[27]=myZ[9]*myZ[27];
   myZ[20]=myZ[22] + myZ[20];
   myZ[20]=myZ[10]*myZ[20];
   myZ[20]=myZ[20] + myZ[23] + myZ[27] + myZ[26];
   myZ[20]=myZ[2]*myZ[20];
   myZ[22]=5./18.*myZ[18] - 4./9. + myZ[24];
   myZ[22]=myZ[8]*myZ[22];
   myZ[23]= - 11./9.*myZ[18] + 20./9. - myZ[19];
   myZ[23]=myZ[5]*myZ[23];
   myZ[24]= - 16./9. + myZ[25];
   myZ[24]=myZ[7]*myZ[24];
   myZ[20]=myZ[24] + myZ[22] + myZ[23] + myZ[20];
   myZ[20]=myZ[6]*myZ[20];
   myZ[22]= - 1./8. - myZ[11];
   myZ[23]=MMH*myZ[11];
   myZ[23]=myZ[13] - myZ[14] + myZ[23];
   myZ[23]=myZ[2]*myZ[23];
   myZ[22]=1./12.*myZ[23] + 1./3.*myZ[22];
   myZ[22]=myZ[22]*myZ[21]*MMH;
   myZ[23]=11./3. - 3*myZ[19];
   myZ[23]=myZ[23]*myZ[19];
   myZ[23]=myZ[23] + 11./3.*myZ[18];
   myZ[23]=myZ[14]*myZ[23];
   myZ[24]=myZ[13]*myZ[21];
   myZ[23]=myZ[23] - myZ[24];
   myZ[24]= - 1 + 3./4.*myZ[19];
   myZ[24]=myZ[24]*myZ[19];
   myZ[25]=5./3.*myZ[18];
   myZ[24]=myZ[25] + 4 + myZ[24];
   myZ[24]=myZ[15]*myZ[24];
   myZ[22]=myZ[24] + 1./4.*myZ[23] + myZ[22];
   myZ[22]=myZ[2]*myZ[22];
   myZ[23]=myZ[13] - myZ[15];
   myZ[24]=3./4.*myZ[17];
   myZ[23]=myZ[24]*myZ[19]*myZ[23];
   myZ[24]=myZ[25] + myZ[19] - 8./3.;
   myZ[25]=4./3.*myZ[4];
   myZ[24]=myZ[24]*myZ[25];
   myZ[25]= - myZ[24] + 5./24.*myZ[18] + 8 - 209./24.*myZ[19];
   myZ[21]=myZ[11]*myZ[21];
   myZ[24]=myZ[5]*myZ[24];
   myZ[18]=1./12.*myZ[18] + 29./3. - 33./4.*myZ[19];
   myZ[18]=myZ[12]*myZ[18];
   myZ[19]=1 + myZ[12];
   myZ[19]=myZ[19]*pow(CW,2);

      return myZ[18] + 4*myZ[19] + myZ[20] + myZ[21] + myZ[22] + 
      myZ[23] + myZ[24] + 1./3.*myZ[25];
}
