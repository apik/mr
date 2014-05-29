#include <ZZ.hpp>
std::complex<long double> ZZ::my11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> myZ[33];

    myZ[1]=pow(CW,-1);
    myZ[2]=pow(MMZ,-1);
    myZ[3]=pow(SW,-1);
    myZ[4]=double(nL + nH);
    myZ[5]=Tsil::B(0,0,MMZ,mu2);
    myZ[6]=prot00000->M(0);
    myZ[7]=Tsil::I2(0,0,MMt,mu2);
    myZ[8]=Tsil::B(MMt,MMt,MMZ,mu2);
    myZ[9]=Tsil::A(MMt,mu2);
    myZ[10]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    myZ[11]=pow(MMt,-1);
    myZ[12]=Tsil::Aeps(MMt,mu2);
    myZ[13]=prottttt0->M(0);
    myZ[14]=prottttt0->Vzxyv(0);
    myZ[15]=prottttt0->Suxv(0);
    myZ[16]=1/(4*MMt - MMZ);
   myZ[17]=pow(myZ[3],2);
   myZ[18]=pow(myZ[1],2);
   myZ[19]=25./9.*myZ[18] + myZ[17] - 64./9.;
   myZ[20]=4*myZ[9];
   myZ[21]=myZ[19]*myZ[20]*myZ[16];
   myZ[22]=25./3.*myZ[18] - 64./3. + 3*myZ[17];
   myZ[23]=myZ[16]*MMZ;
   myZ[22]=myZ[22]*myZ[23];
   myZ[23]=myZ[23] + 1;
   myZ[23]=myZ[19]*myZ[23];
   myZ[24]= - myZ[8]*myZ[23];
   myZ[22]=myZ[24] - myZ[21] + myZ[22] + 143./9.*myZ[18] - 320./9. + 7*
   myZ[17];
   myZ[22]=myZ[8]*myZ[22];
   myZ[24]=17./9.*myZ[18] + myZ[17] - 32./9.;
   myZ[25]=4./3.*myZ[24];
   myZ[26]=myZ[12]*myZ[25];
   myZ[27]= - myZ[15]*myZ[24];
   myZ[26]=myZ[26] + myZ[27];
   myZ[27]=myZ[17] + myZ[18];
   myZ[28]=myZ[11]*pow(myZ[9],2);
   myZ[28]=myZ[28] - myZ[7];
   myZ[28]= - 2*myZ[28];
   myZ[28]=myZ[27]*myZ[28];
   myZ[29]=41./9.*myZ[18] - 128./9. + myZ[17];
   myZ[29]=myZ[9]*myZ[29];
   myZ[30]= - 7./9.*myZ[18] + myZ[17] + 64./9.;
   myZ[31]=myZ[8]*myZ[9];
   myZ[32]=myZ[30]*myZ[31];
   myZ[26]= - 4./3.*myZ[32] + 2*myZ[26] + 7./3.*myZ[29] + myZ[28];
   myZ[28]=29./3.*myZ[18] - 128./3. - myZ[17];
   myZ[29]= - myZ[8]*myZ[30];
   myZ[28]=2*myZ[28] + myZ[29];
   myZ[29]=2*myZ[8];
   myZ[28]=myZ[28]*myZ[29];
   myZ[29]=myZ[13] + 2*myZ[14];
   myZ[29]=MMt*myZ[29];
   myZ[29]=8*myZ[29] + 4*myZ[10];
   myZ[29]=myZ[30]*myZ[29];
   myZ[28]=myZ[28] - 197./18.*myZ[18] - 128./9. - 29./2.*myZ[17] + 
   myZ[29];
   myZ[29]=myZ[31] + myZ[12];
   myZ[20]=myZ[20] - myZ[15] + 2*myZ[29];
   myZ[20]=myZ[27]*myZ[20];
   myZ[29]=2 + myZ[8];
   myZ[29]=myZ[8]*myZ[29];
   myZ[29]=myZ[29] - 1;
   myZ[29]=MMt*myZ[27]*myZ[29];
   myZ[20]=myZ[29] + myZ[20];
   myZ[20]=myZ[2]*myZ[20];
   myZ[20]=4./3.*myZ[20] + 1./3.*myZ[28];
   myZ[20]=MMt*myZ[20];
   myZ[20]=2*myZ[26] + myZ[20];
   myZ[20]=myZ[2]*myZ[20];
   myZ[26]=myZ[6]*MMZ;
   myZ[28]=4./3.*MMZ;
   myZ[28]=myZ[13]*myZ[28];
   myZ[28]= - 4./3.*myZ[26] - 2*myZ[5] + myZ[28];
   myZ[24]=myZ[24]*myZ[28];
   myZ[27]= - myZ[13]*myZ[27];
   myZ[25]= - myZ[14]*myZ[25];
   myZ[25]=myZ[27] + myZ[25];
   myZ[25]=MMt*myZ[25];
   myZ[27]=17./18.*myZ[18] - 16./9. + 1./2.*myZ[17];
   myZ[28]=4*myZ[16];
   myZ[19]= - myZ[12]*myZ[19]*myZ[28];
   myZ[23]= - myZ[10]*myZ[23];
   myZ[17]=11./9.*myZ[18] + myZ[17] - 20./9.;
   myZ[18]=8./3.*myZ[26] + 31./3. + 4*myZ[5];
   myZ[17]=myZ[4]*myZ[17]*myZ[18];

      return myZ[17] + myZ[19] + myZ[20] + myZ[21] + myZ[22] + myZ[23]
       + myZ[24] + 4*myZ[25] + 5*myZ[27];
}
