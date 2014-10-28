#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::m11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> armZZbar[29], mZZbarret;

    armZZbar[1]=double(nH);
    armZZbar[2]=pow(CW,-1);
    armZZbar[3]=pow(MMH,-1);
    armZZbar[4]=pow(MMZ,-1);
    armZZbar[5]=pow(SW,-1);
    armZZbar[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    armZZbar[7]=Tsil::A(MMt,mu2);
    armZZbar[8]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    armZZbar[9]=Tsil::Aeps(MMt,mu2);
    armZZbar[10]=std::real(Tsil::B(0,0,MMZ,mu2));
    armZZbar[11]=prottttt0->M(0);
    armZZbar[12]=prot00000->M(0);
    armZZbar[13]=prottttt0->Vzxyv(0);
    armZZbar[14]=prottttt0->Suxv(0);
    armZZbar[15]=double(nL);
    armZZbar[16]=1/(4*MMt - MMZ);
   armZZbar[17]=2*armZZbar[8];
   armZZbar[18]= - 1 + armZZbar[17];
   armZZbar[19]= - 4 - 1./3.*armZZbar[6];
   armZZbar[19]=armZZbar[6]*armZZbar[19];
   armZZbar[18]=1./3.*armZZbar[18] + armZZbar[19];
   armZZbar[19]= - 2 - armZZbar[6];
   armZZbar[19]=armZZbar[6]*armZZbar[19];
   armZZbar[17]=armZZbar[19] - 35 + armZZbar[17];
   armZZbar[19]=pow(armZZbar[5],2);
   armZZbar[17]=armZZbar[19]*armZZbar[17];
   armZZbar[20]= - 299 - 14*armZZbar[8];
   armZZbar[21]=58 + 7./3.*armZZbar[6];
   armZZbar[21]=armZZbar[6]*armZZbar[21];
   armZZbar[20]=1./3.*armZZbar[20] + armZZbar[21];
   armZZbar[21]=pow(armZZbar[2],2);
   armZZbar[20]=armZZbar[21]*armZZbar[20];
   armZZbar[17]=1./3.*armZZbar[20] + 64./3.*armZZbar[18] + armZZbar[17]
   ;
   armZZbar[18]=4*armZZbar[3];
   armZZbar[20]=armZZbar[18] + 1./3.*armZZbar[13];
   armZZbar[20]=2*armZZbar[20] + 1./3.*armZZbar[11];
   armZZbar[20]=armZZbar[19]*armZZbar[20];
   armZZbar[18]=armZZbar[18] - 7./27.*armZZbar[13];
   armZZbar[18]=2*armZZbar[18] - 7./27.*armZZbar[11];
   armZZbar[18]=armZZbar[21]*armZZbar[18];
   armZZbar[22]=2*armZZbar[13] + armZZbar[11];
   armZZbar[18]=armZZbar[18] + 64./27.*armZZbar[22] + armZZbar[20];
   armZZbar[18]=MMt*armZZbar[18];
   armZZbar[17]=1./3.*armZZbar[17] + 4*armZZbar[18];
   armZZbar[17]=MMt*armZZbar[17];
   armZZbar[18]=2 + armZZbar[6];
   armZZbar[18]=armZZbar[6]*armZZbar[18];
   armZZbar[18]= - 1 + armZZbar[18];
   armZZbar[20]=armZZbar[19]*armZZbar[18];
   armZZbar[18]=armZZbar[21]*armZZbar[18];
   armZZbar[18]=armZZbar[20] + armZZbar[18];
   armZZbar[18]=MMt*armZZbar[18];
   armZZbar[20]=armZZbar[6]*armZZbar[7];
   armZZbar[22]=2*armZZbar[20] + 4*armZZbar[7] - armZZbar[14] + 2*
   armZZbar[9];
   armZZbar[23]=armZZbar[19]*armZZbar[22];
   armZZbar[22]=armZZbar[21]*armZZbar[22];
   armZZbar[18]=armZZbar[18] + armZZbar[23] + armZZbar[22];
   armZZbar[18]=armZZbar[4]*MMt*armZZbar[18];
   armZZbar[22]= - 17*armZZbar[14] + 95./3.*armZZbar[9];
   armZZbar[23]= - 12*armZZbar[7]*armZZbar[3];
   armZZbar[24]=211./27. + armZZbar[23];
   armZZbar[24]=armZZbar[7]*armZZbar[24];
   armZZbar[20]=14./27.*armZZbar[20] + 1./9.*armZZbar[22] + 
   armZZbar[24];
   armZZbar[20]=armZZbar[21]*armZZbar[20];
   armZZbar[22]=11./3. + armZZbar[23];
   armZZbar[22]=armZZbar[7]*armZZbar[22];
   armZZbar[23]= - armZZbar[6]*armZZbar[7];
   armZZbar[22]=2./3.*armZZbar[23] + armZZbar[22] - armZZbar[14] + 7./3.
   *armZZbar[9];
   armZZbar[22]=armZZbar[19]*armZZbar[22];
   armZZbar[23]=4./3.*armZZbar[23] - 14./3.*armZZbar[7] + armZZbar[14]
    - 4./3.*armZZbar[9];
   armZZbar[20]=armZZbar[20] + 32./9.*armZZbar[23] + armZZbar[22];
   armZZbar[17]=2./3.*armZZbar[18] + 2*armZZbar[20] + armZZbar[17];
   armZZbar[17]=armZZbar[4]*armZZbar[17];
   armZZbar[18]= - 4*armZZbar[10];
   armZZbar[20]= - armZZbar[12] - 4*armZZbar[11];
   armZZbar[20]=MMZ*armZZbar[20];
   armZZbar[22]=MMZ*armZZbar[8];
   armZZbar[22]=4*armZZbar[9] + armZZbar[22];
   armZZbar[22]=armZZbar[16]*armZZbar[22];
   armZZbar[23]= - armZZbar[7]*armZZbar[16];
   armZZbar[20]=64*armZZbar[23] + 16*armZZbar[22] + 8./3.*armZZbar[20]
    + 16*armZZbar[8] - 215./3. + armZZbar[18];
   armZZbar[22]=armZZbar[16]*MMZ;
   armZZbar[24]=1 + armZZbar[22];
   armZZbar[24]=armZZbar[6]*armZZbar[24];
   armZZbar[25]=armZZbar[7]*armZZbar[16];
   armZZbar[26]= - armZZbar[16]*MMZ;
   armZZbar[24]=1./3.*armZZbar[24] + 4./3.*armZZbar[25] - 5./3. + 
   armZZbar[26];
   armZZbar[24]=armZZbar[6]*armZZbar[24];
   armZZbar[20]=1./3.*armZZbar[20] + 16*armZZbar[24];
   armZZbar[24]= - 1 + armZZbar[26];
   armZZbar[24]=armZZbar[6]*armZZbar[24];
   armZZbar[26]=armZZbar[24] + 4*armZZbar[23] + 7 + 3*armZZbar[22];
   armZZbar[26]=armZZbar[6]*armZZbar[26];
   armZZbar[27]=armZZbar[12] + armZZbar[11];
   armZZbar[27]=MMZ*armZZbar[27];
   armZZbar[28]= - MMZ*armZZbar[8];
   armZZbar[28]= - 4*armZZbar[9] + armZZbar[28];
   armZZbar[28]=armZZbar[16]*armZZbar[28];
   armZZbar[26]=armZZbar[26] + 4*armZZbar[25] + armZZbar[28] + 4./3.*
   armZZbar[27] - armZZbar[8] + 77./6. + 2*armZZbar[10];
   armZZbar[26]=armZZbar[19]*armZZbar[26];
   armZZbar[27]=5*armZZbar[12] + 17*armZZbar[11];
   armZZbar[27]=MMZ*armZZbar[27];
   armZZbar[25]=100*armZZbar[25] + 25*armZZbar[28] + 4./3.*armZZbar[27]
    - 25*armZZbar[8] + 937./6. + 10*armZZbar[10];
   armZZbar[22]=25./3.*armZZbar[24] + 100./3.*armZZbar[23] + 143./3. + 
   25*armZZbar[22];
   armZZbar[22]=armZZbar[6]*armZZbar[22];
   armZZbar[22]=1./3.*armZZbar[25] + armZZbar[22];
   armZZbar[22]=armZZbar[21]*armZZbar[22];
   armZZbar[23]= - 4./3.*armZZbar[13] - armZZbar[11];
   armZZbar[23]=armZZbar[19]*armZZbar[23];
   armZZbar[24]= - 68./27.*armZZbar[13] - armZZbar[11];
   armZZbar[24]=armZZbar[21]*armZZbar[24];
   armZZbar[23]=armZZbar[24] + 128./27.*armZZbar[13] + armZZbar[23];
   armZZbar[23]=MMt*armZZbar[23];
   armZZbar[17]=2*armZZbar[17] + 4*armZZbar[23] + 1./3.*armZZbar[22] + 
   4./3.*armZZbar[20] + armZZbar[26];
   armZZbar[17]=armZZbar[1]*armZZbar[17];
   armZZbar[18]= - 31./3. + armZZbar[18];
   armZZbar[18]=armZZbar[15]*armZZbar[18];
   armZZbar[20]= - MMZ*armZZbar[15]*armZZbar[12];
   armZZbar[18]=armZZbar[18] + 8./3.*armZZbar[20];
   armZZbar[20]=31./3. + 4*armZZbar[10];
   armZZbar[20]=armZZbar[15]*armZZbar[20];
   armZZbar[22]=MMZ*armZZbar[15]*armZZbar[12];
   armZZbar[20]=armZZbar[20] + 8./3.*armZZbar[22];
   armZZbar[19]=armZZbar[19]*armZZbar[20];
   armZZbar[20]=armZZbar[21]*armZZbar[20];

      mZZbarret = armZZbar[17] + 20./9.*armZZbar[18] + armZZbar[19] + 
      11./9.*armZZbar[20];
      return mZZbarret;
}
