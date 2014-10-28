#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::my11(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuZZ[29], yuZZret;

    aryuZZ[1]=double(nH);
    aryuZZ[2]=pow(CW,-1);
    aryuZZ[3]=pow(MMH,-1);
    aryuZZ[4]=pow(MMZ,-1);
    aryuZZ[5]=pow(SW,-1);
    aryuZZ[6]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuZZ[7]=Tsil::A(MMt,mu2);
    aryuZZ[8]=Tsil::Beps(MMt,MMt,MMZ,mu2);
    aryuZZ[9]=Tsil::Aeps(MMt,mu2);
    aryuZZ[10]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryuZZ[11]=prottttt0->M(0);
    aryuZZ[12]=prot00000->M(0);
    aryuZZ[13]=prottttt0->Vzxyv(0);
    aryuZZ[14]=prottttt0->Suxv(0);
    aryuZZ[15]=double(nL);
    aryuZZ[16]=1/(4*MMt - MMZ);
   aryuZZ[17]=2*aryuZZ[8];
   aryuZZ[18]= - 1 + aryuZZ[17];
   aryuZZ[19]= - 4 - 1./3.*aryuZZ[6];
   aryuZZ[19]=aryuZZ[6]*aryuZZ[19];
   aryuZZ[18]=1./3.*aryuZZ[18] + aryuZZ[19];
   aryuZZ[19]= - 2 - aryuZZ[6];
   aryuZZ[19]=aryuZZ[6]*aryuZZ[19];
   aryuZZ[17]=aryuZZ[19] - 35 + aryuZZ[17];
   aryuZZ[19]=pow(aryuZZ[5],2);
   aryuZZ[17]=aryuZZ[19]*aryuZZ[17];
   aryuZZ[20]= - 299 - 14*aryuZZ[8];
   aryuZZ[21]=58 + 7./3.*aryuZZ[6];
   aryuZZ[21]=aryuZZ[6]*aryuZZ[21];
   aryuZZ[20]=1./3.*aryuZZ[20] + aryuZZ[21];
   aryuZZ[21]=pow(aryuZZ[2],2);
   aryuZZ[20]=aryuZZ[21]*aryuZZ[20];
   aryuZZ[17]=1./3.*aryuZZ[20] + 64./3.*aryuZZ[18] + aryuZZ[17];
   aryuZZ[18]=4*aryuZZ[3];
   aryuZZ[20]=aryuZZ[18] + 1./3.*aryuZZ[13];
   aryuZZ[20]=2*aryuZZ[20] + 1./3.*aryuZZ[11];
   aryuZZ[20]=aryuZZ[19]*aryuZZ[20];
   aryuZZ[18]=aryuZZ[18] - 7./27.*aryuZZ[13];
   aryuZZ[18]=2*aryuZZ[18] - 7./27.*aryuZZ[11];
   aryuZZ[18]=aryuZZ[21]*aryuZZ[18];
   aryuZZ[22]=2*aryuZZ[13] + aryuZZ[11];
   aryuZZ[18]=aryuZZ[18] + 64./27.*aryuZZ[22] + aryuZZ[20];
   aryuZZ[18]=MMt*aryuZZ[18];
   aryuZZ[17]=1./3.*aryuZZ[17] + 4*aryuZZ[18];
   aryuZZ[17]=MMt*aryuZZ[17];
   aryuZZ[18]=2 + aryuZZ[6];
   aryuZZ[18]=aryuZZ[6]*aryuZZ[18];
   aryuZZ[18]= - 1 + aryuZZ[18];
   aryuZZ[20]=aryuZZ[19]*aryuZZ[18];
   aryuZZ[18]=aryuZZ[21]*aryuZZ[18];
   aryuZZ[18]=aryuZZ[20] + aryuZZ[18];
   aryuZZ[18]=MMt*aryuZZ[18];
   aryuZZ[20]=aryuZZ[6]*aryuZZ[7];
   aryuZZ[22]=2*aryuZZ[20] + 4*aryuZZ[7] - aryuZZ[14] + 2*aryuZZ[9];
   aryuZZ[23]=aryuZZ[19]*aryuZZ[22];
   aryuZZ[22]=aryuZZ[21]*aryuZZ[22];
   aryuZZ[18]=aryuZZ[18] + aryuZZ[23] + aryuZZ[22];
   aryuZZ[18]=aryuZZ[4]*MMt*aryuZZ[18];
   aryuZZ[22]= - 17*aryuZZ[14] + 95./3.*aryuZZ[9];
   aryuZZ[23]= - 12*aryuZZ[7]*aryuZZ[3];
   aryuZZ[24]=211./27. + aryuZZ[23];
   aryuZZ[24]=aryuZZ[7]*aryuZZ[24];
   aryuZZ[20]=14./27.*aryuZZ[20] + 1./9.*aryuZZ[22] + aryuZZ[24];
   aryuZZ[20]=aryuZZ[21]*aryuZZ[20];
   aryuZZ[22]=11./3. + aryuZZ[23];
   aryuZZ[22]=aryuZZ[7]*aryuZZ[22];
   aryuZZ[23]= - aryuZZ[6]*aryuZZ[7];
   aryuZZ[22]=2./3.*aryuZZ[23] + aryuZZ[22] - aryuZZ[14] + 7./3.*
   aryuZZ[9];
   aryuZZ[22]=aryuZZ[19]*aryuZZ[22];
   aryuZZ[23]=4./3.*aryuZZ[23] - 14./3.*aryuZZ[7] + aryuZZ[14] - 4./3.*
   aryuZZ[9];
   aryuZZ[20]=aryuZZ[20] + 32./9.*aryuZZ[23] + aryuZZ[22];
   aryuZZ[17]=2./3.*aryuZZ[18] + 2*aryuZZ[20] + aryuZZ[17];
   aryuZZ[17]=aryuZZ[4]*aryuZZ[17];
   aryuZZ[18]= - 4*aryuZZ[10];
   aryuZZ[20]= - aryuZZ[12] - 4*aryuZZ[11];
   aryuZZ[20]=MMZ*aryuZZ[20];
   aryuZZ[22]=MMZ*aryuZZ[8];
   aryuZZ[22]=4*aryuZZ[9] + aryuZZ[22];
   aryuZZ[22]=aryuZZ[16]*aryuZZ[22];
   aryuZZ[23]= - aryuZZ[7]*aryuZZ[16];
   aryuZZ[20]=64*aryuZZ[23] + 16*aryuZZ[22] + 8./3.*aryuZZ[20] + 16*
   aryuZZ[8] - 215./3. + aryuZZ[18];
   aryuZZ[22]=aryuZZ[16]*MMZ;
   aryuZZ[24]=1 + aryuZZ[22];
   aryuZZ[24]=aryuZZ[6]*aryuZZ[24];
   aryuZZ[25]=aryuZZ[7]*aryuZZ[16];
   aryuZZ[26]= - aryuZZ[16]*MMZ;
   aryuZZ[24]=1./3.*aryuZZ[24] + 4./3.*aryuZZ[25] - 5./3. + aryuZZ[26];
   aryuZZ[24]=aryuZZ[6]*aryuZZ[24];
   aryuZZ[20]=1./3.*aryuZZ[20] + 16*aryuZZ[24];
   aryuZZ[24]= - 1 + aryuZZ[26];
   aryuZZ[24]=aryuZZ[6]*aryuZZ[24];
   aryuZZ[26]=aryuZZ[24] + 4*aryuZZ[23] + 7 + 3*aryuZZ[22];
   aryuZZ[26]=aryuZZ[6]*aryuZZ[26];
   aryuZZ[27]=aryuZZ[12] + aryuZZ[11];
   aryuZZ[27]=MMZ*aryuZZ[27];
   aryuZZ[28]= - MMZ*aryuZZ[8];
   aryuZZ[28]= - 4*aryuZZ[9] + aryuZZ[28];
   aryuZZ[28]=aryuZZ[16]*aryuZZ[28];
   aryuZZ[26]=aryuZZ[26] + 4*aryuZZ[25] + aryuZZ[28] + 4./3.*aryuZZ[27]
    - aryuZZ[8] + 77./6. + 2*aryuZZ[10];
   aryuZZ[26]=aryuZZ[19]*aryuZZ[26];
   aryuZZ[27]=5*aryuZZ[12] + 17*aryuZZ[11];
   aryuZZ[27]=MMZ*aryuZZ[27];
   aryuZZ[25]=100*aryuZZ[25] + 25*aryuZZ[28] + 4./3.*aryuZZ[27] - 25*
   aryuZZ[8] + 937./6. + 10*aryuZZ[10];
   aryuZZ[22]=25./3.*aryuZZ[24] + 100./3.*aryuZZ[23] + 143./3. + 25*
   aryuZZ[22];
   aryuZZ[22]=aryuZZ[6]*aryuZZ[22];
   aryuZZ[22]=1./3.*aryuZZ[25] + aryuZZ[22];
   aryuZZ[22]=aryuZZ[21]*aryuZZ[22];
   aryuZZ[23]= - 4./3.*aryuZZ[13] - aryuZZ[11];
   aryuZZ[23]=aryuZZ[19]*aryuZZ[23];
   aryuZZ[24]= - 68./27.*aryuZZ[13] - aryuZZ[11];
   aryuZZ[24]=aryuZZ[21]*aryuZZ[24];
   aryuZZ[23]=aryuZZ[24] + 128./27.*aryuZZ[13] + aryuZZ[23];
   aryuZZ[23]=MMt*aryuZZ[23];
   aryuZZ[17]=2*aryuZZ[17] + 4*aryuZZ[23] + 1./3.*aryuZZ[22] + 4./3.*
   aryuZZ[20] + aryuZZ[26];
   aryuZZ[17]=aryuZZ[1]*aryuZZ[17];
   aryuZZ[18]= - 31./3. + aryuZZ[18];
   aryuZZ[18]=aryuZZ[15]*aryuZZ[18];
   aryuZZ[20]= - MMZ*aryuZZ[15]*aryuZZ[12];
   aryuZZ[18]=aryuZZ[18] + 8./3.*aryuZZ[20];
   aryuZZ[20]=31./3. + 4*aryuZZ[10];
   aryuZZ[20]=aryuZZ[15]*aryuZZ[20];
   aryuZZ[22]=MMZ*aryuZZ[15]*aryuZZ[12];
   aryuZZ[20]=aryuZZ[20] + 8./3.*aryuZZ[22];
   aryuZZ[19]=aryuZZ[19]*aryuZZ[20];
   aryuZZ[20]=aryuZZ[21]*aryuZZ[20];

      yuZZret = aryuZZ[17] + 20./9.*aryuZZ[18] + aryuZZ[19] + 11./9.*
      aryuZZ[20];
      return yuZZret;
}
