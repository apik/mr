#include <ZZ.hpp>
std::complex<long double> zz::m10(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mZZ2MS[27], mZZ2MSret;

    mZZ2MS[1]=pow(c,-1);
    mZZ2MS[2]=pow(mmH,-1);
    mZZ2MS[3]=pow(mmZ,-1);
    mZZ2MS[4]=pow(s,-1);
    mZZ2MS[5]=double(nL + nH);
    mZZ2MS[6]=std::real(Tsil::B(0,0,mmZ,mu2));
    mZZ2MS[7]=double(nH);
    mZZ2MS[8]=Tsil::B(mmt,mmt,mmZ,mu2);
    mZZ2MS[9]=Tsil::A(mmt,mu2);
    mZZ2MS[10]=double(nL);
    mZZ2MS[11]=Tsil::B(mmZ,mmH,mmZ,mu2);
    mZZ2MS[12]=Tsil::B(mmW,mmW,mmZ,mu2);
    mZZ2MS[13]=Tsil::A(mmH,mu2);
    mZZ2MS[14]=Tsil::A(mmZ,mu2);
    mZZ2MS[15]=Tsil::A(mmW,mu2);
   mZZ2MS[16]=pow(mZZ2MS[1],2);
   mZZ2MS[17]=pow(mZZ2MS[4],2);
   mZZ2MS[18]=17./9.*mZZ2MS[16] + mZZ2MS[17] - 32./9.;
   mZZ2MS[19]=1./2.*mZZ2MS[17];
   mZZ2MS[20]= - 7./18.*mZZ2MS[16] + 32./9. + mZZ2MS[19];
   mZZ2MS[20]=mZZ2MS[8]*mZZ2MS[20];
   mZZ2MS[21]=mZZ2MS[16] + mZZ2MS[17];
   mZZ2MS[22]=mZZ2MS[2]*mZZ2MS[9]*mZZ2MS[21];
   mZZ2MS[20]=6*mZZ2MS[22] + mZZ2MS[20] - mZZ2MS[18];
   mZZ2MS[20]=mmt*mZZ2MS[20];
   mZZ2MS[18]= - mZZ2MS[9]*mZZ2MS[18];
   mZZ2MS[18]=mZZ2MS[18] + mZZ2MS[20];
   mZZ2MS[18]=mZZ2MS[7]*mZZ2MS[18];
   mZZ2MS[20]=mZZ2MS[21]*mmH;
   mZZ2MS[22]= - mZZ2MS[13] + mZZ2MS[14];
   mZZ2MS[22]=mZZ2MS[22]*mZZ2MS[20];
   mZZ2MS[23]=mZZ2MS[21]*mZZ2MS[11];
   mZZ2MS[24]= - pow(mmH,2)*mZZ2MS[23];
   mZZ2MS[22]=mZZ2MS[24] + mZZ2MS[22];
   mZZ2MS[22]=mZZ2MS[3]*mZZ2MS[22];
   mZZ2MS[24]= - 1./6.*mZZ2MS[14] - 1./2.*mZZ2MS[13];
   mZZ2MS[24]=mZZ2MS[21]*mZZ2MS[24];
   mZZ2MS[25]=1./3.*mZZ2MS[11] + 1./6.;
   mZZ2MS[20]=mZZ2MS[25]*mZZ2MS[20];
   mZZ2MS[25]=1./6.*mZZ2MS[16];
   mZZ2MS[26]= - mZZ2MS[25] - 4 + 5./2.*mZZ2MS[17];
   mZZ2MS[26]=mZZ2MS[15]*mZZ2MS[26];
   mZZ2MS[18]=1./12.*mZZ2MS[22] + mZZ2MS[18] + mZZ2MS[26] + mZZ2MS[20]
    + mZZ2MS[24];
   mZZ2MS[18]=mZZ2MS[3]*mZZ2MS[18];
   mZZ2MS[20]=11./9.*mZZ2MS[16] + mZZ2MS[17] - 20./9.;
   mZZ2MS[22]= - 17./18.*mZZ2MS[16] + 16./9. - mZZ2MS[19];
   mZZ2MS[22]=mZZ2MS[8]*mZZ2MS[22];
   mZZ2MS[19]= - 5./18.*mZZ2MS[16] + 4./9. - mZZ2MS[19];
   mZZ2MS[19]=mZZ2MS[6]*mZZ2MS[19];
   mZZ2MS[19]=mZZ2MS[19] + 1./3.*mZZ2MS[20] + mZZ2MS[22];
   mZZ2MS[19]=mZZ2MS[7]*mZZ2MS[19];
   mZZ2MS[20]= - mZZ2MS[10]*mZZ2MS[20];
   mZZ2MS[22]=mZZ2MS[17] - 4;
   mZZ2MS[22]=mZZ2MS[16] + 1./3.*mZZ2MS[22];
   mZZ2MS[22]= - mZZ2MS[5]*mZZ2MS[22];
   mZZ2MS[20]=mZZ2MS[22] + mZZ2MS[20];
   mZZ2MS[22]=mZZ2MS[6] - 1./3.;
   mZZ2MS[20]=mZZ2MS[22]*mZZ2MS[20];
   mZZ2MS[22]=3*mZZ2MS[17];
   mZZ2MS[24]= - mZZ2MS[16] + 2 - mZZ2MS[22];
   mZZ2MS[24]=mmZ*mZZ2MS[24];
   mZZ2MS[21]=mZZ2MS[14]*mZZ2MS[21];
   mZZ2MS[22]= - mZZ2MS[15]*mZZ2MS[22];
   mZZ2MS[21]= - 3./2.*mZZ2MS[21] + mZZ2MS[22] + mZZ2MS[24];
   mZZ2MS[21]=mZZ2MS[2]*mZZ2MS[21];
   mZZ2MS[22]=pow(c,2);
   mZZ2MS[22]=4*mZZ2MS[22];
   mZZ2MS[16]= - mZZ2MS[22] - 1./12.*mZZ2MS[16] - 29./3. + 33./4.*
   mZZ2MS[17];
   mZZ2MS[16]=mZZ2MS[12]*mZZ2MS[16];
   mZZ2MS[17]=mZZ2MS[25] - 8 + 59./6.*mZZ2MS[17];

      mZZ2MSret = mZZ2MS[16] + 1./3.*mZZ2MS[17] + mZZ2MS[18] + 
      mZZ2MS[19] + mZZ2MS[20] + mZZ2MS[21] - mZZ2MS[22] - mZZ2MS[23];
      return mZZ2MSret;
}
