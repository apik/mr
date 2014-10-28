#include <ZZ.hpp>
std::complex<long double>
ZZ<OS>::my10(size_t nL, size_t nH, size_t boson)
{     
      
      
    std::complex<long double> aryuZZ[25], yuZZret;

    aryuZZ[1]=double(nL + nH);
    aryuZZ[2]=pow(CW,-1);
    aryuZZ[3]=pow(SW,-1);
    aryuZZ[4]=std::real(Tsil::B(0,0,MMZ,mu2));
    aryuZZ[5]=double(nH);
    aryuZZ[6]=pow(MMZ,-1);
    aryuZZ[7]=Tsil::B(MMt,MMt,MMZ,mu2);
    aryuZZ[8]=Tsil::A(MMt,mu2);
    aryuZZ[9]=pow(MMH,-1);
    aryuZZ[10]=double(nL);
    aryuZZ[11]=double(boson);
    aryuZZ[12]=Tsil::B(MMZ,MMH,MMZ,mu2);
    aryuZZ[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    aryuZZ[14]=Tsil::A(MMH,mu2);
    aryuZZ[15]=Tsil::A(MMZ,mu2);
    aryuZZ[16]=Tsil::A(MMW,mu2);
   aryuZZ[17]=MMH*aryuZZ[12];
   aryuZZ[17]=aryuZZ[17] + aryuZZ[14] - aryuZZ[15];
   aryuZZ[17]=MMH*aryuZZ[17];
   aryuZZ[18]=pow(aryuZZ[2],2);
   aryuZZ[19]=aryuZZ[18]*aryuZZ[17];
   aryuZZ[20]=pow(aryuZZ[3],2);
   aryuZZ[17]=aryuZZ[20]*aryuZZ[17];
   aryuZZ[17]=aryuZZ[19] + aryuZZ[17];
   aryuZZ[17]=aryuZZ[6]*aryuZZ[17];
   aryuZZ[19]=1./3.*aryuZZ[15];
   aryuZZ[21]=aryuZZ[19] + 1./3.*aryuZZ[16] + aryuZZ[14];
   aryuZZ[22]= - 1./2. - aryuZZ[12];
   aryuZZ[22]=1./3.*MMH*aryuZZ[22];
   aryuZZ[21]=1./2.*aryuZZ[21] + aryuZZ[22];
   aryuZZ[21]=aryuZZ[18]*aryuZZ[21];
   aryuZZ[19]=aryuZZ[19] - 5*aryuZZ[16] + aryuZZ[14];
   aryuZZ[19]=1./2.*aryuZZ[19] + aryuZZ[22];
   aryuZZ[19]=aryuZZ[20]*aryuZZ[19];
   aryuZZ[17]=1./12.*aryuZZ[17] + aryuZZ[19] + 4*aryuZZ[16] + 
   aryuZZ[21];
   aryuZZ[17]=aryuZZ[6]*aryuZZ[17];
   aryuZZ[19]= - 1./3. + 1./2.*aryuZZ[13];
   aryuZZ[21]=MMZ + 3./2.*aryuZZ[15];
   aryuZZ[21]=aryuZZ[9]*aryuZZ[21];
   aryuZZ[19]=aryuZZ[21] + 1./6.*aryuZZ[19] + aryuZZ[12];
   aryuZZ[19]=aryuZZ[18]*aryuZZ[19];
   aryuZZ[21]= - 59./9. - 33./2.*aryuZZ[13];
   aryuZZ[22]=1./2.*aryuZZ[15] + MMZ + aryuZZ[16];
   aryuZZ[22]=aryuZZ[9]*aryuZZ[22];
   aryuZZ[21]=3*aryuZZ[22] + 1./2.*aryuZZ[21] + aryuZZ[12];
   aryuZZ[21]=aryuZZ[20]*aryuZZ[21];
   aryuZZ[22]=pow(CW,2);
   aryuZZ[23]=2./3. + aryuZZ[22];
   aryuZZ[22]=29./3. + 4*aryuZZ[22];
   aryuZZ[22]=aryuZZ[13]*aryuZZ[22];
   aryuZZ[24]= - aryuZZ[9]*MMZ;
   aryuZZ[17]=aryuZZ[17] + aryuZZ[21] + aryuZZ[19] + 2*aryuZZ[24] + 4*
   aryuZZ[23] + aryuZZ[22];
   aryuZZ[17]=aryuZZ[11]*aryuZZ[17];
   aryuZZ[19]=5./3.*aryuZZ[10] + aryuZZ[1];
   aryuZZ[21]= - 5./3.*aryuZZ[10] - aryuZZ[1];
   aryuZZ[21]=aryuZZ[4]*aryuZZ[21];
   aryuZZ[22]= - aryuZZ[4] + 5./3. - 4*aryuZZ[7];
   aryuZZ[22]=aryuZZ[5]*aryuZZ[22];
   aryuZZ[19]=1./3.*aryuZZ[22] + 1./3.*aryuZZ[19] + aryuZZ[21];
   aryuZZ[21]= - 11./9.*aryuZZ[10] - aryuZZ[1];
   aryuZZ[22]=11./9.*aryuZZ[10] + aryuZZ[1];
   aryuZZ[22]=aryuZZ[4]*aryuZZ[22];
   aryuZZ[23]=5./2.*aryuZZ[4] - 11./3. + 17./2.*aryuZZ[7];
   aryuZZ[23]=aryuZZ[5]*aryuZZ[23];
   aryuZZ[21]=1./9.*aryuZZ[23] + 1./3.*aryuZZ[21] + aryuZZ[22];
   aryuZZ[21]=aryuZZ[18]*aryuZZ[21];
   aryuZZ[22]= - aryuZZ[10] - 1./3.*aryuZZ[1];
   aryuZZ[23]=aryuZZ[10] + 1./3.*aryuZZ[1];
   aryuZZ[23]=aryuZZ[4]*aryuZZ[23];
   aryuZZ[24]=1./2.*aryuZZ[4] - 1./3. + 1./2.*aryuZZ[7];
   aryuZZ[24]=aryuZZ[5]*aryuZZ[24];
   aryuZZ[22]=aryuZZ[24] + 1./3.*aryuZZ[22] + aryuZZ[23];
   aryuZZ[22]=aryuZZ[20]*aryuZZ[22];
   aryuZZ[23]=17 + 7./2.*aryuZZ[7];
   aryuZZ[23]=MMt*aryuZZ[23];
   aryuZZ[23]=17*aryuZZ[8] + aryuZZ[23];
   aryuZZ[24]= - 6*aryuZZ[9]*MMt*aryuZZ[8];
   aryuZZ[23]=1./9.*aryuZZ[23] + aryuZZ[24];
   aryuZZ[18]=aryuZZ[18]*aryuZZ[5]*aryuZZ[23];
   aryuZZ[23]=1 - 1./2.*aryuZZ[7];
   aryuZZ[23]=MMt*aryuZZ[23];
   aryuZZ[23]=aryuZZ[24] + aryuZZ[8] + aryuZZ[23];
   aryuZZ[20]=aryuZZ[20]*aryuZZ[5]*aryuZZ[23];
   aryuZZ[23]= - 1 - aryuZZ[7];
   aryuZZ[23]=MMt*aryuZZ[23];
   aryuZZ[23]= - aryuZZ[8] + aryuZZ[23];
   aryuZZ[23]=aryuZZ[5]*aryuZZ[23];
   aryuZZ[18]=aryuZZ[20] + 32./9.*aryuZZ[23] + aryuZZ[18];
   aryuZZ[18]=aryuZZ[6]*aryuZZ[18];

      yuZZret = aryuZZ[17] + aryuZZ[18] + 4./3.*aryuZZ[19] + aryuZZ[21]
       + aryuZZ[22];
      return yuZZret;
}
