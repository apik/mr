#include <ZZ.hpp>
std::complex<long double> ZZ::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mZZ[27];

    mZZ[1]=pow(CW,-1);
    mZZ[2]=pow(MMH,-1);
    mZZ[3]=pow(MMZ,-1);
    mZZ[4]=pow(SW,-1);
    mZZ[5]=double(nL + nH);
    mZZ[6]=Tsil::B(0,0,MMZ,mu2);
    mZZ[7]=double(nH);
    mZZ[8]=Tsil::B(MMt,MMt,MMZ,mu2);
    mZZ[9]=Tsil::B(MMb,MMb,MMZ,mu2);
    mZZ[10]=Tsil::A(MMt,mu2);
    mZZ[11]=Tsil::A(MMb,mu2);
    mZZ[12]=Tsil::B(MMZ,MMH,MMZ,mu2);
    mZZ[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    mZZ[14]=Tsil::A(MMH,mu2);
    mZZ[15]=Tsil::A(MMZ,mu2);
    mZZ[16]=Tsil::A(MMW,mu2);
   mZZ[17]=pow(mZZ[1],2);
   mZZ[18]=pow(mZZ[4],2);
   mZZ[19]=mZZ[18] - 8./9. + 5./9.*mZZ[17];
   mZZ[20]=1./2.*mZZ[18];
   mZZ[21]=mZZ[20] + 17./18.*mZZ[17];
   mZZ[22]= - 8./9. - mZZ[21];
   mZZ[22]=mZZ[9]*mZZ[22];
   mZZ[22]=mZZ[22] + mZZ[19];
   mZZ[22]=MMb*mZZ[22];
   mZZ[23]=mZZ[18] - 32./9. + 17./9.*mZZ[17];
   mZZ[24]=7./18.*mZZ[17] - 32./9. - mZZ[20];
   mZZ[24]=mZZ[8]*mZZ[24];
   mZZ[24]=mZZ[24] + mZZ[23];
   mZZ[24]=MMt*mZZ[24];
   mZZ[25]= - mZZ[10]*MMt;
   mZZ[26]= - MMb*mZZ[11];
   mZZ[25]=mZZ[25] + mZZ[26];
   mZZ[26]=mZZ[17] + mZZ[18];
   mZZ[25]=mZZ[2]*mZZ[26]*mZZ[25];
   mZZ[23]=mZZ[10]*mZZ[23];
   mZZ[19]=mZZ[11]*mZZ[19];
   mZZ[19]=6*mZZ[25] + mZZ[19] + mZZ[24] + mZZ[23] + mZZ[22];
   mZZ[19]=mZZ[3]*mZZ[19];
   mZZ[21]= - 16./9. + mZZ[21];
   mZZ[21]=mZZ[8]*mZZ[21];
   mZZ[20]=5./18.*mZZ[17] - 4./9. + mZZ[20];
   mZZ[20]=mZZ[9]*mZZ[20];
   mZZ[22]= - 11./9.*mZZ[17] + 20./9. - mZZ[18];
   mZZ[22]=mZZ[6]*mZZ[22];
   mZZ[19]=mZZ[22] + mZZ[19] + mZZ[21] + mZZ[20];
   mZZ[19]=mZZ[7]*mZZ[19];
   mZZ[20]=mZZ[12]*pow(MMH,2);
   mZZ[21]= - mZZ[15]*MMH;
   mZZ[20]=mZZ[20] + mZZ[21];
   mZZ[20]=mZZ[3]*mZZ[20];
   mZZ[21]= - 1./2. - mZZ[12];
   mZZ[21]=MMH*mZZ[21];
   mZZ[22]=mZZ[3]*MMH;
   mZZ[22]=1./6.*mZZ[22] + 1;
   mZZ[22]=mZZ[14]*mZZ[22];
   mZZ[20]=1./2.*mZZ[22] + 1./12.*mZZ[20] + 1./3.*mZZ[21];
   mZZ[20]=mZZ[26]*mZZ[20];
   mZZ[21]=mZZ[26]*mZZ[15];
   mZZ[22]=1./6.*mZZ[17];
   mZZ[23]=mZZ[22] + 4 - 5./2.*mZZ[18];
   mZZ[23]=mZZ[16]*mZZ[23];
   mZZ[20]=1./6.*mZZ[21] + mZZ[23] + mZZ[20];
   mZZ[20]=mZZ[3]*mZZ[20];
   mZZ[23]=3*mZZ[18];
   mZZ[24]=mZZ[17] - 2 + mZZ[23];
   mZZ[24]=MMZ*mZZ[24];
   mZZ[23]=mZZ[16]*mZZ[23];
   mZZ[21]=mZZ[24] + 3./2.*mZZ[21] + mZZ[23];
   mZZ[21]=mZZ[2]*mZZ[21];
   mZZ[23]=pow(CW,2);
   mZZ[23]=4*mZZ[23];
   mZZ[24]=mZZ[23] + 1./12.*mZZ[17] + 29./3. - 33./4.*mZZ[18];
   mZZ[24]=mZZ[13]*mZZ[24];
   mZZ[22]= - mZZ[22] + 8 - 59./6.*mZZ[18];
   mZZ[25]=mZZ[12]*mZZ[26];
   mZZ[17]=5./3.*mZZ[17] + mZZ[18] - 8./3.;
   mZZ[18]= - 1./3. + mZZ[6];
   mZZ[17]=mZZ[5]*mZZ[17]*mZZ[18];

      return 4./3.*mZZ[17] + mZZ[19] + mZZ[20] + mZZ[21] + 1./3.*
      mZZ[22] + mZZ[23] + mZZ[24] + mZZ[25];
}
