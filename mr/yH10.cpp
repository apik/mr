#include <HH.hpp>
std::complex<long double> HH::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myH[23], myHret;

    myH[1]=pow(CW,-1);
    myH[2]=pow(MMH,-1);
    myH[3]=pow(MMZ,-1);
    myH[4]=pow(SW,-1);
    myH[5]=double(nH);
    myH[6]=Tsil::B(MMt,MMt,MMH,mu2);
    myH[7]=Tsil::A(MMt,mu2);
    myH[8]=Tsil::B(MMH,MMH,MMH,mu2);
    myH[9]=Tsil::B(MMZ,MMZ,MMH,mu2);
    myH[10]=Tsil::B(MMW,MMW,MMH,mu2);
    myH[11]=Tsil::A(MMZ,mu2);
    myH[12]=Tsil::A(MMW,mu2);
    myH[13]=1/( - MMW + MMH);
    myH[14]=Tsil::A(MMH,mu2);
   myH[15]=myH[6] - 1./2.;
   myH[15]=myH[15]*MMt;
   myH[15]=myH[15] - myH[7];
   myH[16]=3./2.*myH[5];
   myH[15]=myH[15]*myH[16];
   myH[16]=9./2.*myH[8] + myH[10] + 1./2.;
   myH[17]=1./4.*MMH;
   myH[16]=myH[16]*myH[17];
   myH[17]=myH[6]*pow(MMt,2)*myH[2]*myH[5];
   myH[18]=myH[9]*MMH;
   myH[15]=myH[15] + myH[16] - 6*myH[17] + 2*myH[12] + myH[11] + 1./8.*
   myH[18];
   myH[16]=pow(myH[4],2);
   myH[17]= - myH[11] + myH[12];
   myH[17]=myH[17]*myH[16];
   myH[17]=3./4.*myH[17] + myH[15];
   myH[17]=myH[17]*myH[16];
   myH[18]=pow(myH[1],2);
   myH[15]=myH[15]*myH[18];
   myH[15]=myH[17] + myH[15];
   myH[15]=myH[3]*myH[15];
   myH[17]=myH[14] - myH[12];
   myH[19]=3./4.*myH[13];
   myH[17]=myH[19]*myH[17];
   myH[19]=myH[2]*MMZ;
   myH[20]=3*myH[19];
   myH[21]=myH[20] - 1;
   myH[22]=1./2.*myH[9];
   myH[21]=myH[21]*myH[22];
   myH[22]= - 1 + myH[10];
   myH[20]=myH[22]*myH[20];
   myH[17]=myH[21] + myH[20] + 3./8. - myH[10] + myH[17];
   myH[16]=myH[16]*myH[17];
   myH[17]=2 - 3*myH[10];
   myH[17]=myH[17]*myH[19];
   myH[19]=myH[21] + 1./8. - myH[19];
   myH[18]=myH[19]*myH[18];

      myHret = myH[15] + myH[16] + myH[17] + myH[18];
      return myHret;
}
