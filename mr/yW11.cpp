#include <WW.hpp>
std::complex<long double> WW::my11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myW[21];

    myW[1]=pow(CW,-1);
    myW[2]=pow(MMZ,-1);
    myW[3]=pow(SW,-1);
    myW[4]=double(nL + nH);
    myW[5]=Tsil::B(0,0,MMW,mu2);
    myW[6]=prot00000->M(0);
    myW[7]=Tsil::I2(0,0,MMt,mu2);
    myW[8]=Tsil::B(0,MMt,MMW,mu2);
    myW[9]=Tsil::A(MMt,mu2);
    myW[10]=pow(MMt,-1);
    myW[11]=Tsil::Aeps(MMt,mu2);
    myW[12]=prot00tt0->M(0);
    myW[13]=prot00tt0->Tuxv(0);
   myW[14]= - 1 + 1./3.*myW[8];
   myW[15]=10*myW[8];
   myW[14]=myW[14]*myW[15];
   myW[14]=8./3.*myW[13] + myW[14] + 23./3.;
   myW[14]=MMt*myW[14];
   myW[15]= - 5 + 2*myW[8];
   myW[16]=2*myW[9];
   myW[17]=myW[16]*myW[10];
   myW[15]=myW[17] + 5./3.*myW[15];
   myW[15]=myW[15]*myW[16];
   myW[17]=4./3.*myW[11];
   myW[14]=4*myW[7] - myW[17] + myW[15] + myW[14];
   myW[15]=pow(myW[1],2);
   myW[18]=pow(myW[3],2);
   myW[19]=myW[15] + myW[18];
   myW[14]=myW[19]*myW[14];
   myW[19]=MMt*myW[12];
   myW[20]= - 2./3.*myW[19] + 3*myW[8];
   myW[20]=myW[20]*MMt;
   myW[17]=myW[20] + myW[17];
   myW[17]=myW[17]*MMt;
   myW[16]= - myW[16] + 5*MMt;
   myW[20]=1./3.*myW[9];
   myW[16]=myW[16]*myW[20];
   myW[20]=4./3.*myW[13];
   myW[20]=myW[20]*pow(MMt,2);
   myW[16]=myW[17] + myW[16] + myW[20];
   myW[17]= - myW[15] - 1;
   myW[15]=myW[15]*myW[17];
   myW[15]=myW[15] - myW[18];
   myW[15]=myW[2]*myW[16]*myW[15];
   myW[14]=2*myW[15] + myW[14];
   myW[14]=myW[2]*myW[14];
   myW[15]=1 + 2./3.*myW[8];
   myW[15]=myW[8]*myW[15];
   myW[16]=myW[4] - 1;
   myW[17]=myW[5]*myW[16];
   myW[15]= - myW[17] + myW[19] - 2./3. - myW[15];
   myW[16]=myW[6]*myW[16];
   myW[16]=myW[16] + myW[12];
   myW[17]=8./3.*MMZ;
   myW[16]=myW[16]*myW[17];
   myW[17]= - 1 + myW[8];
   myW[17]=myW[17]*myW[9];
   myW[17]=myW[11] + myW[17];
   myW[17]=myW[10]*myW[17];
   myW[17]=myW[17] + myW[13];
   myW[15]=myW[16] + 31./3.*myW[4] + 16./3.*myW[17] - 4*myW[15];
   myW[15]=myW[15]*myW[18];

      return myW[14] + myW[15] - myW[16];
}
