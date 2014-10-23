#include <WW.hpp>
std::complex<long double> WW::my11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myW[21], myWret;

    myW[1]=double(nH);
    myW[2]=pow(CW,-1);
    myW[3]=pow(MMZ,-1);
    myW[4]=pow(SW,-1);
    myW[5]=Tsil::I2(0,0,MMt,mu2);
    myW[6]=Tsil::B(0,MMt,MMW,mu2);
    myW[7]=Tsil::A(MMt,mu2);
    myW[8]=pow(MMt,-1);
    myW[9]=Tsil::Aeps(MMt,mu2);
    myW[10]=prot00tt0->M(0);
    myW[11]=prot00tt0->Tuxv(0);
    myW[12]=double(nL);
    myW[13]=std::real(Tsil::B(0,0,MMW,mu2));
    myW[14]=prot00000->M(0);
   myW[15]= - 1 + 1./3.*myW[6];
   myW[15]=MMt*myW[15];
   myW[15]=2./3.*myW[7] + myW[15];
   myW[16]=10*myW[6];
   myW[15]=myW[15]*myW[16];
   myW[16]=pow(myW[7],2);
   myW[17]=myW[16]*myW[8];
   myW[17]=myW[17] + myW[5];
   myW[18]= - 23*MMt + 50*myW[7];
   myW[19]=myW[11]*MMt;
   myW[15]= - myW[15] + 1./3.*myW[18] - 8./3.*myW[19] + 4./3.*myW[9] - 
   4*myW[17];
   myW[17]=pow(myW[4],2);
   myW[18]=pow(myW[2],2);
   myW[19]= - myW[18] - myW[17];
   myW[15]=myW[15]*myW[19];
   myW[19]= - 5*myW[7] - 4*myW[9];
   myW[19]=MMt*myW[19];
   myW[16]=2*myW[16] + myW[19];
   myW[19]=2./3.*myW[10];
   myW[19]=myW[19]*pow(MMt,3);
   myW[20]=3*myW[6] + 4./3.*myW[11];
   myW[20]=myW[20]*pow(MMt,2);
   myW[16]=myW[19] - myW[20] + 1./3.*myW[16];
   myW[19]=myW[18] + 1;
   myW[18]=myW[18]*myW[19]*myW[16];
   myW[16]=myW[16]*myW[17];
   myW[16]=myW[16] + myW[18];
   myW[16]=myW[3]*myW[16];
   myW[15]=myW[15] + 2*myW[16];
   myW[15]=myW[3]*myW[15];
   myW[16]=1 + 2./3.*myW[6];
   myW[16]=myW[6]*myW[16];
   myW[18]= - MMt + 2./3.*MMZ;
   myW[18]=myW[10]*myW[18];
   myW[16]=myW[16] + myW[18];
   myW[18]=myW[6] - 1;
   myW[18]=myW[7]*myW[18];
   myW[18]=myW[9] + myW[18];
   myW[18]=myW[8]*myW[18];
   myW[18]=myW[18] + myW[11];
   myW[16]=13 + 4*myW[16] + 16./3.*myW[18];
   myW[16]=myW[16]*myW[17];
   myW[18]=8./3.*MMZ;
   myW[19]= - myW[10]*myW[18];
   myW[15]=myW[15] + myW[19] + myW[16];
   myW[15]=myW[1]*myW[15];
   myW[16]=myW[18]*myW[14];
   myW[18]=myW[16] + 31./3. + 4*myW[13];
   myW[17]=myW[18]*myW[17];
   myW[16]= - myW[16] + myW[17];
   myW[16]=myW[12]*myW[16];

      myWret = myW[15] + myW[16];
      return myWret;
}
