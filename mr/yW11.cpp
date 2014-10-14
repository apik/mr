#include <WW.hpp>
std::complex<long double> WW::my11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myW[17], myWret;

    myW[1]=pow(CW,-1);
    myW[2]=pow(MMH,-1);
    myW[3]=pow(MMZ,-1);
    myW[4]=pow(SW,-1);
    myW[5]=Tsil::I2(0,0,MMt,mu2);
    myW[6]=Tsil::B(0,MMt,MMW,mu2);
    myW[7]=Tsil::A(MMt,mu2);
    myW[8]=pow(MMt,-1);
    myW[9]=Tsil::Aeps(MMt,mu2);
   myW[10]=3*myW[7];
   myW[11]=myW[10] + MMt;
   myW[11]=myW[11]*myW[7];
   myW[12]=pow(MMt,2);
   myW[11]=myW[11] + 3*myW[12];
   myW[13]=16*myW[2];
   myW[11]=myW[11]*myW[13];
   myW[13]=2*myW[7];
   myW[14]=myW[13]*myW[8];
   myW[14]=myW[14] - 7;
   myW[13]=myW[14]*myW[13];
   myW[14]=myW[10] - MMt;
   myW[15]=4*myW[6];
   myW[15]=myW[14]*myW[15];
   myW[16]=myW[9] - myW[5];
   myW[11]= - myW[11] - myW[13] - myW[15] + 21./2.*MMt - 4*myW[16];
   myW[13]=pow(myW[4],2);
   myW[15]=pow(myW[1],2);
   myW[16]=myW[13] + myW[15];
   myW[11]=myW[11]*myW[16];
   myW[16]= - myW[15] - 1;
   myW[15]=myW[15]*myW[16];
   myW[13]= - myW[13] + myW[15];
   myW[10]=myW[10]*MMt;
   myW[10]=myW[10] - myW[12];
   myW[10]=myW[10]*myW[6];
   myW[12]=myW[14]*myW[7];
   myW[10]=myW[10] + myW[12];
   myW[10]=myW[3]*myW[10]*myW[13];
   myW[10]=4*myW[10] + myW[11];

      myWret = myW[10]*myW[3];
      return myWret;
}
