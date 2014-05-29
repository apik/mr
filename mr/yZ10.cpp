#include <ZZ.hpp>
std::complex<long double> ZZ::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myZ[20];

    myZ[1]=pow(CW,-1);
    myZ[2]=pow(MMZ,-1);
    myZ[3]=pow(SW,-1);
    myZ[4]=double(nL + nH);
    myZ[5]=Tsil::B(0,0,MMZ,mu2);
    myZ[6]=Tsil::B(MMZ,MMH,MMZ,mu2);
    myZ[7]=Tsil::B(MMW,MMW,MMZ,mu2);
    myZ[8]=Tsil::B(MMt,MMt,MMZ,mu2);
    myZ[9]=Tsil::A(MMH,mu2);
    myZ[10]=Tsil::A(MMZ,mu2);
    myZ[11]=Tsil::A(MMW,mu2);
    myZ[12]=Tsil::A(MMt,mu2);
    myZ[13]=1/( - MMW + MMH);
   myZ[14]=MMH*myZ[6];
   myZ[14]=myZ[9] + myZ[14] - myZ[10];
   myZ[15]=myZ[2]*MMH;
   myZ[15]=1./12.*myZ[15];
   myZ[14]=myZ[14]*myZ[15];
   myZ[15]=myZ[6] + 1./8.;
   myZ[16]=1./3.*MMH;
   myZ[15]=myZ[15]*myZ[16];
   myZ[14]= - myZ[14] + myZ[15] - 11./12.*myZ[10] + 1./4.*myZ[9];
   myZ[15]=pow(myZ[3],2);
   myZ[16]=myZ[11] - myZ[10];
   myZ[16]=myZ[16]*myZ[15];
   myZ[17]=myZ[8]*MMt;
   myZ[17]=myZ[17] + myZ[12];
   myZ[16]=3./4.*myZ[16] + 1./4.*MMt - myZ[11] - myZ[14] - 1./2.*
   myZ[17];
   myZ[16]=myZ[2]*myZ[16];
   myZ[18]= - 19./18. - 3*myZ[7];
   myZ[18]=myZ[8] + 11./2.*myZ[18] - myZ[5];
   myZ[19]= - myZ[11] + myZ[9];
   myZ[19]=myZ[13]*myZ[19];
   myZ[16]=3./4.*myZ[19] + myZ[6] + myZ[16] + 1./2.*myZ[18];
   myZ[16]=myZ[15]*myZ[16];
   myZ[18]=41./12.*MMt + 5*myZ[11];
   myZ[14]=1./3.*myZ[18] - myZ[14] + 7./18.*myZ[17];
   myZ[14]=myZ[2]*myZ[14];
   myZ[18]=myZ[5] - 1./3.;
   myZ[18]=myZ[18]*myZ[4];
   myZ[19]=5./6. + myZ[7];
   myZ[19]=1./2.*myZ[19] - 17./3.*myZ[5];
   myZ[14]=20./9.*myZ[18] + myZ[14] + 17./18.*myZ[8] + 1./6.*myZ[19] + 
   myZ[6];
   myZ[14]=myZ[14]*pow(myZ[1],2);
   myZ[19]=pow(CW,2);
   myZ[17]=myZ[17] + MMt;
   myZ[17]=myZ[11] - 8./9.*myZ[17];
   myZ[17]=myZ[2]*myZ[17];
   myZ[17]=myZ[17] + 2./3. + myZ[19];
   myZ[15]= - 8./3. + myZ[15];
   myZ[15]=myZ[15]*myZ[18];
   myZ[18]=myZ[8] - myZ[5];
   myZ[19]=29./3. + 4*myZ[19];
   myZ[19]=myZ[7]*myZ[19];

      return myZ[14] + 4./3.*myZ[15] + myZ[16] + 4*myZ[17] - 16./9.*
      myZ[18] + myZ[19];
}
