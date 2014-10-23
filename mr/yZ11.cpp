#include <ZZ.hpp>
std::complex<long double> ZZ::my11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myZ[29], myZret;

    myZ[1]=double(nH);
    myZ[2]=pow(CW,-1);
    myZ[3]=pow(MMZ,-1);
    myZ[4]=pow(SW,-1);
    myZ[5]=Tsil::I2(0,0,MMt,mu2);
    myZ[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    myZ[7]=Tsil::A(MMt,mu2);
    myZ[8]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    myZ[9]=pow(MMt,-1);
    myZ[10]=Tsil::Aeps(MMt,mu2);
    myZ[11]=std::real(Tsil::B(0,0,MMZ,mu2));
    myZ[12]=prottttt0->M(0);
    myZ[13]=prot00000->M(0);
    myZ[14]=prottttt0->Vzxyv(0);
    myZ[15]=prottttt0->Suxv(0);
    myZ[16]=double(nL);
    myZ[17]=1/(4*MMt - MMZ);
   myZ[18]=pow(myZ[4],2);
   myZ[19]=pow(myZ[2],2);
   myZ[20]= - 7./9.*myZ[19] + myZ[18] + 64./9.;
   myZ[21]=myZ[20]*myZ[6];
   myZ[22]=29./3.*myZ[19] - 128./3. - myZ[18];
   myZ[22]=2*myZ[22] - myZ[21];
   myZ[22]=myZ[6]*myZ[22];
   myZ[23]=4*myZ[20];
   myZ[23]=myZ[8]*myZ[23];
   myZ[24]=myZ[12] + 2*myZ[14];
   myZ[20]=MMt*myZ[20]*myZ[24];
   myZ[20]=8*myZ[20] + 2*myZ[22] + myZ[23] - 197./18.*myZ[19] - 128./9.
    - 29./2.*myZ[18];
   myZ[22]=myZ[18] + myZ[19];
   myZ[23]=myZ[10]*myZ[22];
   myZ[24]=myZ[6] + 2;
   myZ[24]=myZ[24]*myZ[22];
   myZ[25]=myZ[7]*myZ[24];
   myZ[23]=myZ[23] + myZ[25];
   myZ[24]=myZ[6]*myZ[24];
   myZ[24]=myZ[24] - myZ[22];
   myZ[24]=MMt*myZ[24];
   myZ[25]= - myZ[15]*myZ[22];
   myZ[23]=myZ[25] + 2*myZ[23] + myZ[24];
   myZ[23]=myZ[3]*myZ[23];
   myZ[20]=4./3.*myZ[23] + 1./3.*myZ[20];
   myZ[20]=MMt*myZ[20];
   myZ[23]=myZ[5]*myZ[22];
   myZ[24]=17./9.*myZ[19] + myZ[18] - 32./9.;
   myZ[25]=myZ[15]*myZ[24];
   myZ[23]=myZ[23] - myZ[25];
   myZ[25]=41./9.*myZ[19] - 128./9. + myZ[18];
   myZ[21]=7*myZ[25] - 4*myZ[21];
   myZ[25]=myZ[7]*myZ[9]*myZ[22];
   myZ[21]=1./3.*myZ[21] - 2*myZ[25];
   myZ[21]=myZ[7]*myZ[21];
   myZ[25]=myZ[10]*myZ[24];
   myZ[21]=8./3.*myZ[25] + myZ[21];
   myZ[20]=2*myZ[21] + 4*myZ[23] + myZ[20];
   myZ[20]=myZ[3]*myZ[20];
   myZ[21]= - MMZ*myZ[8];
   myZ[23]= - myZ[6] + 1;
   myZ[23]=myZ[7]*myZ[23];
   myZ[21]=4*myZ[23] - 4*myZ[10] + myZ[21];
   myZ[23]=25./9.*myZ[19] + myZ[18] - 64./9.;
   myZ[21]=myZ[23]*myZ[21];
   myZ[25]=myZ[23]*myZ[6];
   myZ[26]= - myZ[25] + 25./3.*myZ[19] - 64./3. + 3*myZ[18];
   myZ[26]=myZ[6]*MMZ*myZ[26];
   myZ[21]=myZ[26] + myZ[21];
   myZ[21]=myZ[17]*myZ[21];
   myZ[25]= - myZ[25] + 143./9.*myZ[19] - 320./9. + 7*myZ[18];
   myZ[25]=myZ[6]*myZ[25];
   myZ[22]= - myZ[12]*myZ[22];
   myZ[24]=4./3.*myZ[24];
   myZ[26]= - myZ[14]*myZ[24];
   myZ[22]=myZ[22] + myZ[26];
   myZ[22]=MMt*myZ[22];
   myZ[26]=937./18.*myZ[19] - 860./9. + 77./2.*myZ[18];
   myZ[23]= - myZ[8]*myZ[23];
   myZ[24]=MMZ*myZ[12]*myZ[24];
   myZ[27]=myZ[18] - 8./9. + 5./9.*myZ[19];
   myZ[28]=myZ[11]*myZ[27];
   myZ[20]=myZ[20] + 4*myZ[22] + 2*myZ[28] + myZ[21] + myZ[25] + 
   myZ[24] + 1./3.*myZ[26] + myZ[23];
   myZ[20]=myZ[1]*myZ[20];
   myZ[18]=11./9.*myZ[19] + myZ[18] - 20./9.;
   myZ[18]=myZ[18]*myZ[16];
   myZ[19]=myZ[1]*myZ[27];
   myZ[19]=2*myZ[18] + myZ[19];
   myZ[19]=myZ[13]*MMZ*myZ[19];
   myZ[21]=31./3. + 4*myZ[11];
   myZ[18]=myZ[21]*myZ[18];

      myZret = myZ[18] + 4./3.*myZ[19] + myZ[20];
      return myZret;
}
