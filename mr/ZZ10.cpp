#include <ZZ.hpp>
std::complex<long double> ZZ::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mZZ[26], mZZret;

    mZZ[1]=pow(CW,-1);
    mZZ[2]=pow(MMH,-1);
    mZZ[3]=pow(MMZ,-1);
    mZZ[4]=pow(SW,-1);
    mZZ[5]=double(nL + nH);
    mZZ[6]=std::real(Tsil::B(0,0,MMZ,mu2));
    mZZ[7]=Tsil::B(MMZ,MMH,MMZ,mu2);
    mZZ[8]=Tsil::B(MMW,MMW,MMZ,mu2);
    mZZ[9]=Tsil::B(MMt,MMt,MMZ,mu2);
    mZZ[10]=Tsil::A(MMH,mu2);
    mZZ[11]=Tsil::A(MMZ,mu2);
    mZZ[12]=Tsil::A(MMW,mu2);
    mZZ[13]=Tsil::A(MMt,mu2);
   mZZ[14]= - mZZ[11] + mZZ[10];
   mZZ[14]=MMH*mZZ[14];
   mZZ[15]=mZZ[7]*pow(MMH,2);
   mZZ[14]=mZZ[15] + mZZ[14];
   mZZ[14]=mZZ[3]*mZZ[14];
   mZZ[15]= - 1./3.*mZZ[7] - 1./6.;
   mZZ[15]=MMH*mZZ[15];
   mZZ[14]=1./12.*mZZ[14] + 1./6.*mZZ[11] + 1./2.*mZZ[10] + mZZ[15];
   mZZ[15]=pow(mZZ[1],2);
   mZZ[16]=pow(mZZ[4],2);
   mZZ[17]=mZZ[15] + mZZ[16];
   mZZ[14]=mZZ[17]*mZZ[14];
   mZZ[18]=17./9.*mZZ[15] + mZZ[16] - 32./9.;
   mZZ[19]=mZZ[13]*mZZ[18];
   mZZ[20]=1./6.*mZZ[15];
   mZZ[21]=mZZ[20] + 4 - 5./2.*mZZ[16];
   mZZ[21]=mZZ[12]*mZZ[21];
   mZZ[14]=mZZ[21] + mZZ[19] + mZZ[14];
   mZZ[14]=mZZ[3]*mZZ[14];
   mZZ[19]=1./2.*mZZ[16];
   mZZ[21]=mZZ[19] - 16./9. + 17./18.*mZZ[15];
   mZZ[19]=7./18.*mZZ[15] - 32./9. - mZZ[19];
   mZZ[22]=MMt*mZZ[3];
   mZZ[19]=mZZ[19]*mZZ[22];
   mZZ[19]=mZZ[19] + mZZ[21];
   mZZ[19]=mZZ[9]*mZZ[19];
   mZZ[23]=mZZ[13]*mZZ[22];
   mZZ[23]= - 6*mZZ[23] + 3./2.*mZZ[11];
   mZZ[23]=mZZ[17]*mZZ[23];
   mZZ[24]=3*mZZ[16];
   mZZ[25]=mZZ[15] - 2 + mZZ[24];
   mZZ[25]=MMZ*mZZ[25];
   mZZ[24]=mZZ[12]*mZZ[24];
   mZZ[23]=mZZ[24] + mZZ[25] + mZZ[23];
   mZZ[23]=mZZ[2]*mZZ[23];
   mZZ[20]= - mZZ[20] + 8 - 59./6.*mZZ[16];
   mZZ[17]=mZZ[7]*mZZ[17];
   mZZ[21]= - mZZ[6]*mZZ[21];
   mZZ[18]=mZZ[18]*mZZ[22];
   mZZ[22]=1./12.*mZZ[15] + 29./3. - 33./4.*mZZ[16];
   mZZ[22]=mZZ[8]*mZZ[22];
   mZZ[15]=5./3.*mZZ[15] + mZZ[16] - 8./3.;
   mZZ[16]= - 1./3. + mZZ[6];
   mZZ[15]=mZZ[5]*mZZ[15]*mZZ[16];
   mZZ[16]=1 + mZZ[8];
   mZZ[16]=mZZ[16]*pow(CW,2);

      mZZret = mZZ[14] + 4./3.*mZZ[15] + 4*mZZ[16] + mZZ[17] + mZZ[18]
       + mZZ[19] + 1./3.*mZZ[20] + mZZ[21] + mZZ[22] + mZZ[23];
      return mZZret;
}
