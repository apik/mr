#include <WW.hpp>
std::complex<long double> WW::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myW[28];

    myW[1]=pow(CW,-1);
    myW[2]=pow(MMZ,-1);
    myW[3]=pow(SW,-1);
    myW[4]=double(nL + nH);
    myW[5]=Tsil::B(0,0,MMW,mu2);
    myW[6]=double(nH);
    myW[7]=Tsil::B(MMt,MMb,MMW,mu2);
    myW[8]=Tsil::A(MMt,mu2);
    myW[9]=Tsil::A(MMb,mu2);
    myW[10]=Tsil::B(MMW,MMH,MMW,mu2);
    myW[11]=Tsil::B(MMW,MMZ,MMW,mu2);
    myW[12]=Tsil::A(MMH,mu2);
    myW[13]=Tsil::A(MMZ,mu2);
    myW[14]=Tsil::A(MMW,mu2);
    myW[15]=1/( - MMb + MMt);
    myW[16]=1/( - MMW + MMH);
   myW[17]=1./2.*MMb;
   myW[18]= - myW[9] - myW[17] + 1./2.*MMt;
   myW[18]=myW[18]*myW[6];
   myW[19]=myW[8]*myW[6];
   myW[18]=myW[18] + myW[19];
   myW[20]=3./2.*myW[15];
   myW[18]=myW[20]*MMb*myW[18];
   myW[20]=MMb + MMt;
   myW[21]=myW[7]*myW[6];
   myW[22]=1./2.*myW[21];
   myW[20]=myW[20]*myW[22];
   myW[22]=myW[10] + 1./8.;
   myW[23]=1./3.*MMH;
   myW[22]=myW[22]*myW[23];
   myW[23]=myW[9] + MMb + 1./4.*MMt;
   myW[23]=myW[23]*myW[6];
   myW[18]=myW[18] + myW[20] - myW[23] + myW[22] + 1./2.*myW[19];
   myW[20]=myW[13] - myW[14];
   myW[22]=pow(myW[1],2);
   myW[23]=1./12.*myW[22];
   myW[24]=myW[20]*myW[23];
   myW[25]=53./3.*myW[14] - myW[12];
   myW[24]=myW[24] + 3./2.*myW[13] + 1./4.*myW[25] - myW[18];
   myW[24]=myW[24]*myW[22];
   myW[25]=pow(myW[3],2);
   myW[20]=myW[20]*myW[25];
   myW[26]=5./3.*myW[14] - myW[12];
   myW[18]= - 3./4.*myW[20] - 1./2.*myW[13] + 1./4.*myW[26] - myW[18];
   myW[18]=myW[18]*myW[25];
   myW[20]= - myW[6]*myW[9];
   myW[19]=myW[19] + myW[20];
   myW[20]=MMb - MMt;
   myW[19]=myW[20]*myW[19];
   myW[20]=myW[12] - myW[14];
   myW[26]=MMH*myW[10];
   myW[26]=myW[26] + myW[20];
   myW[27]=1./6.*MMH;
   myW[26]=myW[26]*myW[27];
   myW[19]=myW[26] + myW[19];
   myW[17]=myW[17] - MMt;
   myW[17]=myW[17]*MMb;
   myW[26]=pow(MMt,2);
   myW[17]=myW[17] + 1./2.*myW[26];
   myW[17]=myW[17]*myW[21];
   myW[17]= - myW[17] + 1./2.*myW[19];
   myW[19]=myW[22] + 1;
   myW[19]=myW[22]*myW[19];
   myW[19]=myW[19] + myW[25];
   myW[17]=myW[2]*myW[17]*myW[19];
   myW[17]=myW[17] + myW[24] + myW[18];
   myW[17]=myW[2]*myW[17];
   myW[18]=3./4.*myW[20];
   myW[18]=myW[16]*myW[18];
   myW[19]=4./3.*myW[4] - myW[6];
   myW[19]=myW[5]*myW[19];
   myW[18]=myW[19] + myW[21] - 33./4.*myW[11] - 4./9.*myW[4] - 209./72.
    + myW[10] + myW[18];
   myW[18]=myW[25]*myW[18];
   myW[19]=myW[22] + 17;
   myW[19]=myW[11]*myW[19];
   myW[19]= - 1./2. + myW[19];
   myW[19]=myW[19]*myW[23];
   myW[20]= - 1 + myW[11];

      return myW[17] + myW[18] + myW[19] + 4*myW[20];
}
