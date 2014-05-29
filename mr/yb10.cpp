#include <bb.hpp>
std::complex<long double> bb::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myb[24];

    myb[1]=pow(CW,-1);
    myb[2]=pow(MMZ,-1);
    myb[3]=pow(SW,-1);
    myb[4]=double(nH);
    myb[5]=Tsil::A(MMt,mu2);
    myb[6]=Tsil::B(MMH,MMb,MMb,mu2);
    myb[7]=Tsil::B(MMZ,MMb,MMb,mu2);
    myb[8]=pow(MMb,-1);
    myb[9]=Tsil::B(MMW,MMt,MMb,mu2);
    myb[10]=Tsil::A(MMH,mu2);
    myb[11]=Tsil::A(MMZ,mu2);
    myb[12]=Tsil::A(MMW,mu2);
    myb[13]=Tsil::A(MMb,mu2);
    myb[14]=1/( - MMb + MMt);
    myb[15]=1/( - MMW + MMH);
   myb[16]=MMZ*myb[7];
   myb[16]=myb[13] - myb[11] - myb[16];
   myb[17]=1./2.*myb[8];
   myb[16]=myb[17]*myb[16];
   myb[17]=1./2.*MMt;
   myb[18]=myb[17]*myb[8];
   myb[19]=MMZ*myb[8];
   myb[20]=myb[18] + 1./2. - myb[19];
   myb[20]=myb[9]*myb[20];
   myb[21]=myb[5] - myb[12];
   myb[21]=myb[8]*myb[21];
   myb[22]= - 3./2. + myb[7];
   myb[23]= - myb[12] + myb[10];
   myb[23]=myb[15]*myb[23];
   myb[16]=3./2.*myb[23] + myb[20] + 1./2.*myb[22] + myb[16] + myb[21];
   myb[20]=myb[17] + myb[5];
   myb[22]=1./2.*MMb;
   myb[23]= - myb[13] + myb[20] - myb[22];
   myb[23]=myb[23]*MMb*myb[14];
   myb[20]= - myb[23] - myb[20];
   myb[23]=3./2.*myb[4];
   myb[20]=myb[23]*myb[20];
   myb[17]=myb[21]*myb[17];
   myb[21]=1./2. - myb[6];
   myb[21]=MMH*myb[21];
   myb[21]=myb[5] - myb[10] + 3*myb[11] + 7*myb[12] + myb[21];
   myb[17]=myb[17] + myb[13] + 1./2.*myb[21];
   myb[18]=myb[18] - 1;
   myb[18]=myb[18]*MMt;
   myb[18]=myb[18] + myb[22];
   myb[21]=1./2.*myb[9];
   myb[18]=myb[18]*myb[21];
   myb[21]=MMb*myb[6];
   myb[17]=myb[18] + myb[21] + 1./2.*myb[17] + myb[20];
   myb[18]=pow(myb[3],2);
   myb[20]=myb[12] - myb[11];
   myb[20]=myb[20]*myb[18];
   myb[20]=3./4.*myb[20] + myb[17];
   myb[20]=myb[2]*myb[20];
   myb[16]=1./2.*myb[16] + myb[20];
   myb[16]=myb[16]*myb[18];
   myb[17]=myb[2]*myb[17];
   myb[18]=myb[19]*myb[7];
   myb[20]=myb[11]*myb[8];
   myb[18]=myb[18] + myb[20];
   myb[20]=myb[13]*myb[8];
   myb[21]=myb[18] - myb[20];
   myb[22]= - 1./2. + myb[7];
   myb[21]=17*myb[22] - 5*myb[21];
   myb[17]=1./36.*myb[21] + myb[17];
   myb[17]=myb[17]*pow(myb[1],2);
   myb[16]=myb[16] + myb[17];
   myb[17]=myb[20] - 1 + myb[7];
   myb[17]=myb[18] + 2*myb[17];
   myb[18]=myb[9]*myb[19];

      return 1./2.*myb[16] + 1./9.*myb[17] + 1./4.*myb[18];
}
