#include <bb.hpp>
std::complex<long double> bb::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myb[18];

    myb[1]=pow(CW,-1);
    myb[2]=pow(MMZ,-1);
    myb[3]=pow(SW,-1);
    myb[4]=Tsil::A(MMZ,mu2);
    myb[5]=Tsil::A(MMW,mu2);
    myb[6]=Tsil::A(MMt,mu2);
    myb[7]=Tsil::A(MMb,mu2);
    myb[8]=pow(MMb,-1);
    myb[9]=1/(MMt - MMW);
    myb[10]=1/( - MMW + MMH);
    myb[11]=Tsil::A(MMH,mu2);
   myb[12]=pow(myb[1],2);
   myb[13]= - 3*myb[5] - 1./4.*MMH + 5./4.*MMt + 3./2.*myb[6];
   myb[14]=5./6.*myb[4] - myb[13];
   myb[14]=myb[14]*myb[12];
   myb[15]=pow(myb[3],2);
   myb[16]=myb[5] - myb[4];
   myb[16]=myb[16]*myb[15];
   myb[16]=myb[16] + myb[4];
   myb[13]= - myb[13] + 3./2.*myb[16];
   myb[13]=myb[13]*myb[15];
   myb[13]=myb[14] + myb[13];
   myb[13]= - 1./6. - 1./3.*myb[4] + 1./4.*myb[13];
   myb[13]=myb[2]*myb[13];
   myb[14]= - myb[5] + myb[11];
   myb[14]=myb[10]*myb[14];
   myb[16]=myb[5] - myb[6];
   myb[16]=myb[9]*myb[16];
   myb[16]=myb[16] + 1;
   myb[16]=MMZ*myb[16];
   myb[17]=myb[5] + myb[16];
   myb[17]=myb[9]*myb[17];
   myb[14]=myb[14] + myb[17];
   myb[14]=myb[14]*myb[15];
   myb[15]=myb[9]*myb[16];
   myb[14]=myb[15] - myb[14];
   myb[15]=myb[8]*myb[7];

      return  - 11./72.*myb[12] + myb[13] - 3./8.*myb[14] + 1./3.*
      myb[15];
}
