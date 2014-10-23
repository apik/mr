#include <ZZ.hpp>
std::complex<long double> ZZ::m10(size_t nL, size_t nH, int boson)
{     
      
      
    std::complex<long double> mZZ[26], mZZret;

    mZZ[1]=double(nL + nH);
    mZZ[2]=pow(CW,-1);
    mZZ[3]=pow(SW,-1);
    mZZ[4]=std::real(Tsil::B(0,0,MMZ,mu2));
    mZZ[5]=double(nH);
    mZZ[6]=pow(MMZ,-1);
    mZZ[7]=Tsil::B(MMt,MMt,MMZ,mu2);
    mZZ[8]=Tsil::A(MMt,mu2);
    mZZ[9]=pow(MMH,-1);
    mZZ[10]=double(nL);
    mZZ[11]=double(boson);
    mZZ[12]=Tsil::B(MMZ,MMH,MMZ,mu2);
    mZZ[13]=Tsil::B(MMW,MMW,MMZ,mu2);
    mZZ[14]=Tsil::A(MMH,mu2);
    mZZ[15]=Tsil::A(MMZ,mu2);
    mZZ[16]=Tsil::A(MMW,mu2);
   mZZ[17]= - 1./2. - mZZ[12];
   mZZ[18]=MMH*mZZ[12];
   mZZ[18]=mZZ[18] - mZZ[15] + mZZ[14];
   mZZ[18]=mZZ[6]*mZZ[18];
   mZZ[17]=1./12.*mZZ[18] + 1./3.*mZZ[17];
   mZZ[18]=pow(mZZ[2],2);
   mZZ[19]=pow(mZZ[3],2);
   mZZ[20]=mZZ[18] + mZZ[19];
   mZZ[17]=MMH*mZZ[20]*mZZ[17];
   mZZ[21]=1./2.*mZZ[14] + 1./6.*mZZ[15];
   mZZ[21]=mZZ[20]*mZZ[21];
   mZZ[22]=1./6.*mZZ[18];
   mZZ[23]=mZZ[22] + 4 - 5./2.*mZZ[19];
   mZZ[23]=mZZ[16]*mZZ[23];
   mZZ[17]=mZZ[17] + mZZ[23] + mZZ[21];
   mZZ[17]=mZZ[6]*mZZ[17];
   mZZ[21]=mZZ[15]*mZZ[20];
   mZZ[23]=3*mZZ[19];
   mZZ[24]=mZZ[18] - 2 + mZZ[23];
   mZZ[24]=MMZ*mZZ[24];
   mZZ[23]=mZZ[16]*mZZ[23];
   mZZ[21]=mZZ[24] + 3./2.*mZZ[21] + mZZ[23];
   mZZ[21]=mZZ[9]*mZZ[21];
   mZZ[23]= - 1./3. + 1./2.*mZZ[13];
   mZZ[22]=mZZ[23]*mZZ[22];
   mZZ[20]=mZZ[12]*mZZ[20];
   mZZ[23]=8 - 59./6.*mZZ[19];
   mZZ[24]=29./3. - 33./4.*mZZ[19];
   mZZ[24]=mZZ[13]*mZZ[24];
   mZZ[25]=1 + mZZ[13];
   mZZ[25]=mZZ[25]*pow(CW,2);
   mZZ[17]=4*mZZ[25] + mZZ[17] + mZZ[20] + mZZ[22] + mZZ[24] + 1./3.*
   mZZ[23] + mZZ[21];
   mZZ[17]=mZZ[11]*mZZ[17];
   mZZ[20]= - 32./9. + mZZ[19];
   mZZ[21]=mZZ[8] + MMt;
   mZZ[20]=mZZ[21]*mZZ[20];
   mZZ[22]=6*MMt;
   mZZ[22]=mZZ[22]*mZZ[9]*mZZ[8];
   mZZ[23]= - mZZ[19]*mZZ[22];
   mZZ[21]=17./9.*mZZ[21] - mZZ[22];
   mZZ[21]=mZZ[21]*mZZ[18];
   mZZ[22]=1./2.*mZZ[19];
   mZZ[24]=7./18.*mZZ[18] - 32./9. - mZZ[22];
   mZZ[24]=mZZ[7]*MMt*mZZ[24];
   mZZ[20]=mZZ[24] + mZZ[21] + mZZ[23] + mZZ[20];
   mZZ[20]=mZZ[6]*mZZ[20];
   mZZ[21]=11./9.*mZZ[18] + mZZ[19] - 20./9.;
   mZZ[23]=5./18.*mZZ[18] - 4./9. + mZZ[22];
   mZZ[23]=mZZ[4]*mZZ[23];
   mZZ[22]=17./18.*mZZ[18] - 16./9. + mZZ[22];
   mZZ[22]=mZZ[7]*mZZ[22];
   mZZ[20]=mZZ[20] + mZZ[22] - 1./3.*mZZ[21] + mZZ[23];
   mZZ[20]=mZZ[5]*mZZ[20];
   mZZ[21]=mZZ[10]*mZZ[21];
   mZZ[19]=mZZ[19] - 4;
   mZZ[18]=mZZ[18] + 1./3.*mZZ[19];
   mZZ[18]=mZZ[1]*mZZ[18];
   mZZ[18]=mZZ[18] + mZZ[21];
   mZZ[19]=mZZ[4] - 1./3.;
   mZZ[18]=mZZ[19]*mZZ[18];

      mZZret = mZZ[17] + mZZ[18] + mZZ[20];
      return mZZret;
}
