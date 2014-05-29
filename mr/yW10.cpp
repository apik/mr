#include <WW.hpp>
std::complex<long double> WW::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myW[24];

    myW[1]=pow(CW,-1);
    myW[2]=pow(MMZ,-1);
    myW[3]=pow(SW,-1);
    myW[4]=double(nL + nH);
    myW[5]=Tsil::B(0,0,MMW,mu2);
    myW[6]=Tsil::B(MMW,MMH,MMW,mu2);
    myW[7]=Tsil::B(MMW,MMZ,MMW,mu2);
    myW[8]=Tsil::B(0,MMt,MMW,mu2);
    myW[9]=Tsil::A(MMH,mu2);
    myW[10]=Tsil::A(MMZ,mu2);
    myW[11]=Tsil::A(MMW,mu2);
    myW[12]=Tsil::A(MMt,mu2);
    myW[13]=1/( - MMW + MMH);
   myW[14]=MMH*myW[6];
   myW[15]=myW[9] - myW[11];
   myW[14]=myW[14] + myW[15];
   myW[16]=1./6.*MMH;
   myW[14]=myW[14]*myW[16];
   myW[16]=myW[8]*pow(MMt,2);
   myW[17]=myW[12]*MMt;
   myW[14]= - myW[14] + myW[16] + myW[17];
   myW[16]=1./2.*myW[2];
   myW[14]=myW[14]*myW[16];
   myW[16]=pow(myW[1],2);
   myW[17]= - myW[16] - 1;
   myW[17]=myW[17]*myW[14];
   myW[18]=myW[6] + 1./8.;
   myW[19]=1./3.*MMH;
   myW[18]=myW[18]*myW[19];
   myW[18]=myW[18] + 1./2.*myW[12];
   myW[19]=myW[11] - myW[10];
   myW[20]=1./12.*myW[16];
   myW[21]= - myW[19]*myW[20];
   myW[22]=myW[8] - 1./2.;
   myW[22]=myW[22]*MMt;
   myW[22]=myW[22] + 1./2.*myW[9];
   myW[23]=3*myW[10] + 53./6.*myW[11] - myW[22];
   myW[17]=myW[17] + myW[21] + 1./2.*myW[23] - myW[18];
   myW[17]=myW[2]*myW[16]*myW[17];
   myW[21]=pow(myW[3],2);
   myW[19]=myW[19]*myW[21];
   myW[22]= - myW[10] + 5./6.*myW[11] - myW[22];
   myW[14]=3./4.*myW[19] - myW[14] + 1./2.*myW[22] - myW[18];
   myW[14]=myW[2]*myW[14];
   myW[18]= - 209./8. - 4*myW[4];
   myW[15]=myW[13]*myW[15];
   myW[19]= - 1 + 4./3.*myW[4];
   myW[19]=myW[5]*myW[19];
   myW[14]=myW[19] - 33./4.*myW[7] + 3./4.*myW[15] + myW[6] + 1./9.*
   myW[18] + myW[8] + myW[14];
   myW[14]=myW[14]*myW[21];
   myW[15]=myW[16] + 17;
   myW[15]=myW[7]*myW[15];
   myW[15]= - 1./2. + myW[15];
   myW[15]=myW[15]*myW[20];
   myW[16]= - 1 + myW[7];

      return myW[14] + myW[15] + 4*myW[16] + myW[17];
}
