#include <WW.hpp>
std::complex<long double> ww::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mWW2MS[25], mWW2MSret;

    mWW2MS[1]=pow(c,-1);
    mWW2MS[2]=pow(mmH,-1);
    mWW2MS[3]=pow(mmZ,-1);
    mWW2MS[4]=pow(s,-1);
    mWW2MS[5]=double(nL + nH);
    mWW2MS[6]=std::real(Tsil::B(0,0,mmW,mu2));
    mWW2MS[7]=double(nH);
    mWW2MS[8]=Tsil::B(0,mmt,mmW,mu2);
    mWW2MS[9]=Tsil::A(mmt,mu2);
    mWW2MS[10]=double(nL);
    mWW2MS[11]=Tsil::B(mmW,mmH,mmW,mu2);
    mWW2MS[12]=Tsil::B(mmW,mmZ,mmW,mu2);
    mWW2MS[13]=Tsil::A(mmH,mu2);
    mWW2MS[14]=Tsil::A(mmZ,mu2);
    mWW2MS[15]=Tsil::A(mmW,mu2);
   mWW2MS[16]=mWW2MS[11]*pow(mmH,2);
   mWW2MS[17]=mWW2MS[13]*mmH;
   mWW2MS[18]=mWW2MS[15]*mmH;
   mWW2MS[16]= - mWW2MS[18] + mWW2MS[16] + mWW2MS[17];
   mWW2MS[17]=mWW2MS[8]*mWW2MS[7];
   mWW2MS[18]=mWW2MS[17]*mmt;
   mWW2MS[19]=mWW2MS[9]*mWW2MS[7];
   mWW2MS[18]=mWW2MS[18] + mWW2MS[19];
   mWW2MS[18]=mWW2MS[18]*mmt;
   mWW2MS[16]= - mWW2MS[18] + 1./6.*mWW2MS[16];
   mWW2MS[18]=1./2.*mWW2MS[3];
   mWW2MS[18]=mWW2MS[16]*mWW2MS[18];
   mWW2MS[20]=mWW2MS[19]*mWW2MS[2];
   mWW2MS[20]=6*mWW2MS[20] - mWW2MS[7] + 1./2.*mWW2MS[17];
   mWW2MS[20]=mWW2MS[20]*mmt;
   mWW2MS[21]=1./3.*mmH;
   mWW2MS[22]=mWW2MS[21] - mWW2MS[13];
   mWW2MS[21]=mWW2MS[21]*mWW2MS[11];
   mWW2MS[18]=mWW2MS[18] - mWW2MS[20] - mWW2MS[21] + mWW2MS[19] - 1./2.
   *mWW2MS[22];
   mWW2MS[19]=13./12.*mWW2MS[15] + 5./4.*mWW2MS[14] - mWW2MS[18];
   mWW2MS[19]=mWW2MS[3]*mWW2MS[19];
   mWW2MS[20]=mWW2MS[2]*mmZ;
   mWW2MS[21]=mWW2MS[15]*mWW2MS[2];
   mWW2MS[21]=mWW2MS[20] + mWW2MS[21];
   mWW2MS[22]=mWW2MS[6] - 1./3.;
   mWW2MS[23]=mWW2MS[5]*mWW2MS[22];
   mWW2MS[23]= - mWW2MS[23] + 59./6. + mWW2MS[7];
   mWW2MS[24]=mWW2MS[14]*mWW2MS[2];
   mWW2MS[24]=3./2.*mWW2MS[24];
   mWW2MS[22]= - mWW2MS[10]*mWW2MS[22];
   mWW2MS[17]=mWW2MS[19] - mWW2MS[11] - mWW2MS[24] + mWW2MS[22] + 33./4.
   *mWW2MS[12] - mWW2MS[17] - 3*mWW2MS[21] + 1./3.*mWW2MS[23];
   mWW2MS[17]=mWW2MS[17]*pow(mWW2MS[4],2);
   mWW2MS[18]= - 35./12.*mWW2MS[15] - 3./4.*mWW2MS[14] - mWW2MS[18];
   mWW2MS[18]=mWW2MS[3]*mWW2MS[18];
   mWW2MS[16]= - mWW2MS[3]*mWW2MS[16];
   mWW2MS[19]= - mWW2MS[14] + mWW2MS[15];
   mWW2MS[16]=1./6.*mWW2MS[19] + mWW2MS[16];
   mWW2MS[16]=mWW2MS[3]*mWW2MS[16];
   mWW2MS[16]= - 1./6.*mWW2MS[12] + mWW2MS[16];
   mWW2MS[19]=pow(mWW2MS[1],2);
   mWW2MS[16]=mWW2MS[16]*mWW2MS[19];
   mWW2MS[21]=1 - 17./2.*mWW2MS[12];
   mWW2MS[16]=1./2.*mWW2MS[16] + mWW2MS[18] - mWW2MS[24] + 1./6.*
   mWW2MS[21] - mWW2MS[20];
   mWW2MS[16]=mWW2MS[16]*mWW2MS[19];
   mWW2MS[18]=1 - mWW2MS[12];
   mWW2MS[18]=2*mWW2MS[18] + mWW2MS[20];

      mWW2MSret = mWW2MS[16] + mWW2MS[17] + 2*mWW2MS[18];
      return mWW2MSret;
}
