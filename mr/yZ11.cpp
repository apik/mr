#include <ZZ.hpp>
std::complex<long double> ZZ::my11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myZ[20], myZret;

    myZ[1]=pow(CW,-1);
    myZ[2]=pow(MMH,-1);
    myZ[3]=pow(MMZ,-1);
    myZ[4]=pow(SW,-1);
    myZ[5]=Tsil::I2(0,0,MMt,mu2);
    myZ[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    myZ[7]=Tsil::A(MMt,mu2);
    myZ[8]=pow(MMt,-1);
    myZ[9]=Tsil::Aeps(MMt,mu2);
    myZ[10]=1/(4*MMt - MMZ);
   myZ[11]=pow(myZ[4],2);
   myZ[12]=pow(myZ[1],2);
   myZ[13]=myZ[11] + myZ[12];
   myZ[14]=myZ[13]*myZ[2];
   myZ[15]=3*myZ[11];
   myZ[16]=myZ[15] - 64./3. + 25./3.*myZ[12];
   myZ[17]=myZ[10]*myZ[16];
   myZ[17]= - 6*myZ[14] + myZ[17];
   myZ[18]= - myZ[8]*myZ[13];
   myZ[17]=2*myZ[17] + myZ[18];
   myZ[18]=2*myZ[7];
   myZ[17]=myZ[17]*myZ[18];
   myZ[14]=MMt*myZ[14];
   myZ[19]=7./3.*myZ[12] - 64./3. - myZ[15];
   myZ[19]=myZ[6]*myZ[19];
   myZ[15]=myZ[17] + myZ[19] - 8*myZ[14] + 59./9.*myZ[12] - 128./9. + 
   myZ[15];
   myZ[15]=myZ[15]*myZ[18];
   myZ[17]= - 7./9.*myZ[12] + 64./9. + myZ[11];
   myZ[17]=myZ[6]*myZ[17];
   myZ[14]=2*myZ[17] - 48*myZ[14] + 161./18.*myZ[12] + 128./9. + 25./2.
   *myZ[11];
   myZ[14]=MMt*myZ[14];
   myZ[17]=myZ[5] - myZ[9];
   myZ[13]=myZ[13]*myZ[17];
   myZ[13]=myZ[15] + 4*myZ[13] + myZ[14];
   myZ[13]=myZ[3]*myZ[13];
   myZ[14]=1./2.*myZ[11] - 32./9. + 25./18.*myZ[12];
   myZ[15]= - myZ[6] + 1;
   myZ[17]=myZ[10]*MMZ;
   myZ[17]=myZ[17] + 1;
   myZ[14]=myZ[15]*myZ[17]*myZ[14];
   myZ[11]= - 25./9.*myZ[12] + 64./9. - myZ[11];
   myZ[12]=myZ[6]*myZ[16];
   myZ[11]=4*myZ[11] + myZ[12];
   myZ[11]=myZ[18]*myZ[10]*myZ[11];

      myZret = myZ[11] + myZ[13] + myZ[14];
      return myZret;
}
