#include <HH.hpp>
std::complex<long double> HH::my11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myH[20], myHret;

    myH[1]=double(nH);
    myH[2]=pow(CW,-1);
    myH[3]=pow(MMH,-1);
    myH[4]=pow(MMZ,-1);
    myH[5]=pow(SW,-1);
    myH[6]=Tsil::I2(0,0,MMt,mu2);
    myH[7]=Tsil::B(MMt,MMt,MMH,mu2);
    myH[8]=Tsil::A(MMt,mu2);
    myH[9]=Tsil::Beps(MMt,MMt,MMH,mu2);
    myH[10]=pow(MMt,-1);
    myH[11]=Tsil::Aeps(MMt,mu2);
    myH[12]=prottttt0->M(0);
    myH[13]=prottttt0->Vxzuv(0);
    myH[14]=prottttt0->Suxv(0);
   myH[15]=myH[12] + 2*myH[13];
   myH[16]=8*MMt;
   myH[15]=myH[16]*myH[3]*myH[15];
   myH[16]=10 + 7*myH[7];
   myH[16]=myH[16]*myH[7];
   myH[16]=myH[16] + 25 - 4*myH[9];
   myH[16]=myH[16]*myH[3];
   myH[15]=myH[16] - myH[15] + 4*myH[13] + 6*myH[12];
   myH[16]=4*MMt;
   myH[15]=myH[15]*myH[16];
   myH[17]=myH[8]*myH[7];
   myH[17]= - 18*myH[17] - myH[14] + 6*myH[11];
   myH[18]=4*myH[3];
   myH[17]=myH[17]*myH[18];
   myH[18]=4 + 3*myH[7];
   myH[19]=2*myH[7];
   myH[18]=myH[18]*myH[19];
   myH[15]=myH[15] - myH[17] - myH[18] - 85./2.;
   myH[15]=myH[15]*MMt;
   myH[17]=myH[10]*pow(myH[8],2);
   myH[17]=myH[11] + myH[17] - myH[6];
   myH[18]= - 5 + 6*myH[7];
   myH[18]=myH[18]*myH[8];
   myH[17]= - myH[18] + 2*myH[17];
   myH[16]=myH[12]*myH[16]*MMH;
   myH[15]= - myH[16] + myH[15] + 2*myH[17];
   myH[16]=pow(myH[2],2);
   myH[17]= - pow(myH[5],2);
   myH[16]=myH[17] - myH[16];

      myHret = myH[16]*myH[15]*myH[4]*myH[1];
      return myHret;
}
