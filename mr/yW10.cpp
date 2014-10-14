#include <WW.hpp>
std::complex<long double> WW::my10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myW[23], myWret;

    myW[1]=pow(CW,-1);
    myW[2]=pow(MMZ,-1);
    myW[3]=pow(SW,-1);
    myW[4]=double(nL + nH);
    myW[5]=std::real(Tsil::B(0,0,MMW,mu2));
    myW[6]=Tsil::B(MMW,MMH,MMW,mu2);
    myW[7]=Tsil::B(MMW,MMZ,MMW,mu2);
    myW[8]=Tsil::B(0,MMt,MMW,mu2);
    myW[9]=Tsil::A(MMH,mu2);
    myW[10]=Tsil::A(MMZ,mu2);
    myW[11]=Tsil::A(MMW,mu2);
    myW[12]=Tsil::A(MMt,mu2);
    myW[13]=1/( - MMW + MMH);
   myW[14]=myW[6] + 1./8.;
   myW[15]=1./3.*MMH;
   myW[14]=myW[14]*myW[15];
   myW[15]=MMt - myW[9];
   myW[16]=myW[8]*MMt;
   myW[15]=myW[16] + myW[12] - 1./2.*myW[15];
   myW[16]= - myW[10] - myW[15];
   myW[17]=myW[11] - myW[10];
   myW[18]=pow(myW[3],2);
   myW[19]=myW[17]*myW[18];
   myW[16]=3./4.*myW[19] + 5./12.*myW[11] + 1./2.*myW[16] - myW[14];
   myW[16]=myW[16]*myW[18];
   myW[15]=3*myW[10] - myW[15];
   myW[19]=pow(myW[1],2);
   myW[20]=1./12.*myW[19];
   myW[17]= - myW[17]*myW[20];
   myW[14]=myW[17] + 53./12.*myW[11] + 1./2.*myW[15] - myW[14];
   myW[14]=myW[14]*myW[19];
   myW[15]=myW[11] - myW[9];
   myW[17]=MMH*myW[6];
   myW[17]=myW[17] - myW[15];
   myW[21]=1./6.*MMH;
   myW[17]=myW[17]*myW[21];
   myW[21]=myW[8]*pow(MMt,2);
   myW[22]=MMt*myW[12];
   myW[17]=myW[17] - myW[21] - myW[22];
   myW[21]=myW[19] + 1;
   myW[21]=myW[19]*myW[21];
   myW[21]=myW[18] + myW[21];
   myW[17]=myW[2]*myW[17]*myW[21];
   myW[14]=1./2.*myW[17] + myW[16] + myW[14];
   myW[14]=myW[2]*myW[14];
   myW[16]=3./4.*myW[13];
   myW[15]= - myW[16]*myW[15];
   myW[16]= - 1 + 4./3.*myW[4];
   myW[16]=myW[5]*myW[16];
   myW[15]=myW[16] + myW[8] - 33./4.*myW[7] - 4./9.*myW[4] - 209./72.
    + myW[6] + myW[15];
   myW[15]=myW[18]*myW[15];
   myW[16]=myW[19] + 17;
   myW[16]=myW[7]*myW[16];
   myW[16]= - 1./2. + myW[16];
   myW[16]=myW[16]*myW[20];
   myW[17]= - 1 + myW[7];

      myWret = myW[14] + myW[15] + myW[16] + 4*myW[17];
      return myWret;
}
