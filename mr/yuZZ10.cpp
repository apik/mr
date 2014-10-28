#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuZZ[26], yuZZret;

    aryuZZ[1]=double(nL + nH);
    aryuZZ[2]=pow(CW,-1);
    aryuZZ[3]=pow(SW,-1);
    aryuZZ[4]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryuZZ[5]=double(nH);
    aryuZZ[6]=pow(MMZ,-1);
    aryuZZ[7]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuZZ[8]=Tsil::A(MMt,mu2);
    aryuZZ[9]=double(nL);
    aryuZZ[10]=double(boson);
    aryuZZ[11]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryuZZ[12]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryuZZ[13]=Tsil::A(MMH,mu2);
    aryuZZ[14]=Tsil::A(MMZ,mu2);
    aryuZZ[15]=Tsil::A(MMW,mu2);
    aryuZZ[16]=1/( - MMW + MMH);
   aryuZZ[17]=pow(aryuZZ[3],2);
   aryuZZ[18]=11./3. - 3*aryuZZ[17];
   aryuZZ[18]=aryuZZ[18]*aryuZZ[17];
   aryuZZ[19]=pow(aryuZZ[2],2);
   aryuZZ[18]=aryuZZ[18] + 11./3.*aryuZZ[19];
   aryuZZ[18]=aryuZZ[14]*aryuZZ[18];
   aryuZZ[20]=aryuZZ[19] + aryuZZ[17];
   aryuZZ[21]=aryuZZ[13]*aryuZZ[20];
   aryuZZ[18]=aryuZZ[18] - aryuZZ[21];
   aryuZZ[21]=aryuZZ[20]*MMH;
   aryuZZ[22]=aryuZZ[13] - aryuZZ[14];
   aryuZZ[22]=aryuZZ[22]*aryuZZ[21];
   aryuZZ[20]=pow(MMH,2)*aryuZZ[11]*aryuZZ[20];
   aryuZZ[20]=aryuZZ[20] + aryuZZ[22];
   aryuZZ[20]=aryuZZ[6]*aryuZZ[20];
   aryuZZ[22]= - 1 + 3./4.*aryuZZ[17];
   aryuZZ[22]=aryuZZ[22]*aryuZZ[17];
   aryuZZ[22]=5./3.*aryuZZ[19] + 4 + aryuZZ[22];
   aryuZZ[22]=aryuZZ[15]*aryuZZ[22];
   aryuZZ[23]=aryuZZ[11] + 1./8.;
   aryuZZ[21]= - aryuZZ[23]*aryuZZ[21];
   aryuZZ[18]=1./12.*aryuZZ[20] + aryuZZ[22] + 1./3.*aryuZZ[21] + 1./4.
   *aryuZZ[18];
   aryuZZ[18]=aryuZZ[10]*aryuZZ[18];
   aryuZZ[20]=aryuZZ[8]*aryuZZ[5];
   aryuZZ[21]=aryuZZ[7]*aryuZZ[5];
   aryuZZ[22]=MMt*aryuZZ[21];
   aryuZZ[20]=aryuZZ[22] + aryuZZ[20];
   aryuZZ[22]=1./2.*aryuZZ[17];
   aryuZZ[23]= - aryuZZ[22] - 32./9. + 7./18.*aryuZZ[19];
   aryuZZ[20]=aryuZZ[23]*aryuZZ[20];
   aryuZZ[23]=41./36.*aryuZZ[19] - 32./9. + 1./4.*aryuZZ[17];
   aryuZZ[23]=aryuZZ[23]*MMt*aryuZZ[5];
   aryuZZ[18]=aryuZZ[23] + aryuZZ[18] + aryuZZ[20];
   aryuZZ[18]=aryuZZ[6]*aryuZZ[18];
   aryuZZ[20]=pow(CW,2);
   aryuZZ[20]=4*aryuZZ[20];
   aryuZZ[23]=aryuZZ[20] + 1./12.*aryuZZ[19] + 29./3. - 33./4.*
   aryuZZ[17];
   aryuZZ[23]=aryuZZ[12]*aryuZZ[23];
   aryuZZ[24]= - aryuZZ[15] + aryuZZ[13];
   aryuZZ[24]=aryuZZ[16]*aryuZZ[24];
   aryuZZ[24]=3./4.*aryuZZ[24] - 209./72. + aryuZZ[11];
   aryuZZ[24]=aryuZZ[17]*aryuZZ[24];
   aryuZZ[25]=5./72. + aryuZZ[11];
   aryuZZ[25]=aryuZZ[25]*aryuZZ[19];
   aryuZZ[20]=aryuZZ[23] + aryuZZ[20] + aryuZZ[25] + 8./3. + aryuZZ[24]
   ;
   aryuZZ[20]=aryuZZ[10]*aryuZZ[20];
   aryuZZ[23]=5./18.*aryuZZ[19] - 4./9. + aryuZZ[22];
   aryuZZ[23]=aryuZZ[5]*aryuZZ[23];
   aryuZZ[24]=11./9.*aryuZZ[19] + aryuZZ[17] - 20./9.;
   aryuZZ[25]=aryuZZ[9]*aryuZZ[24];
   aryuZZ[23]=aryuZZ[23] + aryuZZ[25];
   aryuZZ[23]=aryuZZ[4]*aryuZZ[23];
   aryuZZ[25]= - aryuZZ[5] - aryuZZ[9];
   aryuZZ[24]=aryuZZ[24]*aryuZZ[25];
   aryuZZ[22]=17./18.*aryuZZ[19] - 16./9. + aryuZZ[22];
   aryuZZ[21]=aryuZZ[22]*aryuZZ[21];
   aryuZZ[17]=aryuZZ[17] - 4;
   aryuZZ[17]=aryuZZ[19] + 1./3.*aryuZZ[17];
   aryuZZ[19]= - 1./3. + aryuZZ[4];
   aryuZZ[17]=aryuZZ[1]*aryuZZ[17]*aryuZZ[19];

      yuZZret = aryuZZ[17] + aryuZZ[18] + aryuZZ[20] + aryuZZ[21] + 
      aryuZZ[23] + 1./3.*aryuZZ[24];
      return yuZZret;
}
