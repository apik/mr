#include <WW.hpp>
std::complex<long double> ww::m10(size_t nL, size_t nH, int boson)
{     
      
      
    std::complex<long double> mWW2MS[28], mWW2MSret;

    mWW2MS[1]=double(nL + nH);
    mWW2MS[2]=pow(s,-1);
    mWW2MS[3]=std::real(Tsil::B(0,0,mmW,mu2));
    mWW2MS[4]=double(nH);
    mWW2MS[5]=pow(c,-1);
    mWW2MS[6]=pow(mmZ,-1);
    mWW2MS[7]=Tsil::B(0,mmt,mmW,mu2);
    mWW2MS[8]=Tsil::A(mmt,mu2);
    mWW2MS[9]=pow(mmH,-1);
    mWW2MS[10]=double(nL);
    mWW2MS[11]=double(boson);
    mWW2MS[12]=Tsil::B(mmW,mmH,mmW,mu2);
    mWW2MS[13]=Tsil::B(mmW,mmZ,mmW,mu2);
    mWW2MS[14]=Tsil::A(mmH,mu2);
    mWW2MS[15]=Tsil::A(mmZ,mu2);
    mWW2MS[16]=Tsil::A(mmW,mu2);
   mWW2MS[17]=pow(mWW2MS[5],2);
   mWW2MS[18]=mWW2MS[17]*mWW2MS[6];
   mWW2MS[19]= - 3 - 1./3.*mWW2MS[17];
   mWW2MS[19]=mWW2MS[19]*mWW2MS[18];
   mWW2MS[20]=pow(mWW2MS[2],2);
   mWW2MS[21]=mWW2MS[20]*mWW2MS[6];
   mWW2MS[19]=mWW2MS[19] + 5*mWW2MS[21];
   mWW2MS[22]=mWW2MS[17] + mWW2MS[20];
   mWW2MS[23]=3*mWW2MS[9];
   mWW2MS[24]= - mWW2MS[22]*mWW2MS[23];
   mWW2MS[19]=1./2.*mWW2MS[19] + mWW2MS[24];
   mWW2MS[19]=mWW2MS[15]*mWW2MS[19];
   mWW2MS[22]=mWW2MS[22]*mWW2MS[6];
   mWW2MS[24]=mWW2MS[17] + 1;
   mWW2MS[24]=mWW2MS[17]*mWW2MS[24];
   mWW2MS[24]=mWW2MS[24] + mWW2MS[20];
   mWW2MS[24]=mWW2MS[24]*pow(mWW2MS[6],2);
   mWW2MS[25]=1./6.*mmH;
   mWW2MS[26]= - mWW2MS[24]*mWW2MS[25];
   mWW2MS[26]=mWW2MS[26] - mWW2MS[22];
   mWW2MS[26]=mWW2MS[14]*mWW2MS[26];
   mWW2MS[19]=mWW2MS[19] + mWW2MS[26];
   mWW2MS[26]=mWW2MS[24]*mmH;
   mWW2MS[27]= - 35 + mWW2MS[17];
   mWW2MS[18]=mWW2MS[27]*mWW2MS[18];
   mWW2MS[18]=mWW2MS[26] + mWW2MS[18] + 13*mWW2MS[21];
   mWW2MS[21]= - mWW2MS[20]*mWW2MS[23];
   mWW2MS[18]=1./12.*mWW2MS[18] + mWW2MS[21];
   mWW2MS[18]=mWW2MS[16]*mWW2MS[18];
   mWW2MS[21]= - 1./4.*mWW2MS[26] + mWW2MS[22];
   mWW2MS[21]=mmH*mWW2MS[21];
   mWW2MS[21]= - mWW2MS[20] + 1./3.*mWW2MS[21];
   mWW2MS[21]=mWW2MS[12]*mWW2MS[21];
   mWW2MS[23]=mWW2MS[22]*mWW2MS[25];
   mWW2MS[25]= - 17 - mWW2MS[17];
   mWW2MS[25]=mWW2MS[25]*mWW2MS[17];
   mWW2MS[25]=33./4.*mWW2MS[20] - 4 + 1./12.*mWW2MS[25];
   mWW2MS[25]=mWW2MS[13]*mWW2MS[25];
   mWW2MS[26]= - 3*mWW2MS[20] + 2 - mWW2MS[17];
   mWW2MS[26]=mmZ*mWW2MS[9]*mWW2MS[26];
   mWW2MS[17]=mWW2MS[26] + mWW2MS[25] + mWW2MS[21] + mWW2MS[18] + 
   mWW2MS[23] + 59./18.*mWW2MS[20] + 4 + 1./6.*mWW2MS[17] + 1./2.*
   mWW2MS[19];
   mWW2MS[17]=mWW2MS[11]*mWW2MS[17];
   mWW2MS[18]=mmt*mWW2MS[24];
   mWW2MS[18]=mWW2MS[18] + mWW2MS[22];
   mWW2MS[18]=mmt*mWW2MS[18];
   mWW2MS[18]= - mWW2MS[20] + 1./2.*mWW2MS[18];
   mWW2MS[18]=mWW2MS[7]*mWW2MS[18];
   mWW2MS[19]=mWW2MS[9]*mWW2MS[22];
   mWW2MS[19]=1./2.*mWW2MS[24] + 6*mWW2MS[19];
   mWW2MS[19]=mmt*mWW2MS[19];
   mWW2MS[19]=mWW2MS[19] - mWW2MS[22];
   mWW2MS[19]=mWW2MS[8]*mWW2MS[19];
   mWW2MS[21]= - mmt*mWW2MS[22];
   mWW2MS[18]=mWW2MS[19] + mWW2MS[18] + 1./3.*mWW2MS[20] + mWW2MS[21];
   mWW2MS[18]=mWW2MS[4]*mWW2MS[18];
   mWW2MS[19]= - 1./3.*mWW2MS[3] + 1./9.;
   mWW2MS[19]=mWW2MS[1]*mWW2MS[19];
   mWW2MS[21]=1./3. - mWW2MS[3];
   mWW2MS[21]=mWW2MS[10]*mWW2MS[21];
   mWW2MS[19]=mWW2MS[21] + mWW2MS[19];
   mWW2MS[19]=mWW2MS[20]*mWW2MS[19];

      mWW2MSret = mWW2MS[17] + mWW2MS[18] + mWW2MS[19];
      return mWW2MSret;
}
