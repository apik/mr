#include <WW.hpp>
std::complex<long double> WW::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myW[26], myWret;

    myW[1]=pow(CW,-1);
    myW[2]=pow(MMZ,-1);
    myW[3]=pow(SW,-1);
    myW[4]=double(nL + nH);
    myW[5]=std::real(Tsil::B(0,0,MMW,mu2));
    myW[6]=double(nH);
    myW[7]=Tsil::B(0,MMt,MMW,mu2);
    myW[8]=Tsil::A(MMt,mu2);
    myW[9]=double(nL);
    myW[10]=Tsil::B(MMW,MMH,MMW,mu2);
    myW[11]=Tsil::B(MMW,MMZ,MMW,mu2);
    myW[12]=Tsil::A(MMH,mu2);
    myW[13]=Tsil::A(MMZ,mu2);
    myW[14]=Tsil::A(MMW,mu2);
    myW[15]=1/( - MMW + MMH);
   myW[16]=myW[14] - myW[12];
   myW[17]=MMH*myW[10];
   myW[17]=myW[17] - myW[16];
   myW[17]=MMH*myW[17];
   myW[18]=MMt*myW[6];
   myW[19]=myW[18]*myW[7];
   myW[20]=myW[8]*myW[6];
   myW[19]=myW[19] + myW[20];
   myW[19]=myW[19]*MMt;
   myW[17]= - myW[19] + 1./6.*myW[17];
   myW[19]=1./2.*myW[2];
   myW[21]=myW[17]*myW[19];
   myW[22]=myW[7] - 1./2.;
   myW[18]=1./2.*myW[18];
   myW[18]=myW[22]*myW[18];
   myW[22]=myW[10] + 1./8.;
   myW[23]=1./3.*MMH;
   myW[22]=myW[22]*myW[23];
   myW[18]=myW[21] - 1./2.*myW[20] - myW[18] - myW[22];
   myW[20]=myW[14] - myW[13];
   myW[21]=pow(myW[3],2);
   myW[22]=myW[20]*myW[21];
   myW[23]=1./2.*myW[12];
   myW[24]= - myW[13] - myW[23];
   myW[22]=3./4.*myW[22] + 5./12.*myW[14] + 1./2.*myW[24] + myW[18];
   myW[22]=myW[2]*myW[22];
   myW[24]=3./4.*myW[15];
   myW[16]= - myW[24]*myW[16];
   myW[24]= - 1./3. + myW[7];
   myW[24]=myW[6]*myW[24];
   myW[25]=1./3.*myW[4] + myW[9];
   myW[25]=myW[5]*myW[25];
   myW[16]=myW[22] + myW[25] + myW[24] - 1./3.*myW[9] - 1./9.*myW[4] - 
   209./72. + myW[10] + myW[16];
   myW[16]=myW[21]*myW[16];
   myW[22]=3*myW[13] - myW[23];
   myW[18]=53./12.*myW[14] + 1./2.*myW[22] + myW[18];
   myW[18]=myW[2]*myW[18];
   myW[17]=myW[2]*myW[17];
   myW[17]= - 1./6.*myW[20] + myW[17];
   myW[20]=pow(myW[1],2);
   myW[17]=myW[17]*myW[20]*myW[19];
   myW[17]=myW[17] - 1./24. + myW[18];
   myW[17]=myW[17]*myW[20];
   myW[18]=17 + myW[20];
   myW[18]=myW[18]*myW[20];
   myW[18]=1./12.*myW[18] + 4 - 33./4.*myW[21];
   myW[18]=myW[11]*myW[18];

      myWret =  - 4 + myW[16] + myW[17] + myW[18];
      return myWret;
}
