#include <WW.hpp>
std::complex<long double> ww::m11(size_t nL, size_t nH, int boson)
{     
      
      
    std::complex<long double> mWW2MS[23], mWW2MSret;

    mWW2MS[1]=double(nH);
    mWW2MS[2]=pow(c,-1);
    mWW2MS[3]=pow(mmH,-1);
    mWW2MS[4]=pow(mmZ,-1);
    mWW2MS[5]=pow(s,-1);
    mWW2MS[6]=Tsil::B(0,mmt,mmW,mu2);
    mWW2MS[7]=Tsil::A(mmt,mu2);
    mWW2MS[8]=pow(mmt,-1);
    mWW2MS[9]=Tsil::Aeps(mmt,mu2);
    mWW2MS[10]=prot00tt0->M(0);
    mWW2MS[11]=prot00tt0->Tuxv(0);
    mWW2MS[12]=double(nL);
    mWW2MS[13]=std::real(Tsil::B(0,0,mmW,mu2));
    mWW2MS[14]=prot00000->M(0);
   mWW2MS[15]= - 7 + 5./3.*mWW2MS[6];
   mWW2MS[16]=2*mWW2MS[6];
   mWW2MS[15]=mWW2MS[15]*mWW2MS[16];
   mWW2MS[15]=mWW2MS[15] - 17./6.;
   mWW2MS[15]=mWW2MS[15]*mmt;
   mWW2MS[16]=mWW2MS[7] + 3*mmt;
   mWW2MS[16]=mWW2MS[16]*mmt;
   mWW2MS[17]=pow(mWW2MS[7],2);
   mWW2MS[16]=mWW2MS[16] + 3*mWW2MS[17];
   mWW2MS[16]= - mWW2MS[3]*mWW2MS[16];
   mWW2MS[18]=mWW2MS[11]*mmt;
   mWW2MS[18]=mWW2MS[18] + mWW2MS[9];
   mWW2MS[19]=mWW2MS[6]*mWW2MS[7];
   mWW2MS[20]= - 14*mWW2MS[19] + 23*mWW2MS[7];
   mWW2MS[21]=mWW2MS[17]*mWW2MS[8];
   mWW2MS[15]=16*mWW2MS[16] - mWW2MS[15] - 8*mWW2MS[21] + 4./3.*
   mWW2MS[20] - 8./3.*mWW2MS[18];
   mWW2MS[16]=pow(mWW2MS[2],2);
   mWW2MS[18]=pow(mWW2MS[5],2);
   mWW2MS[20]=mWW2MS[16] + mWW2MS[18];
   mWW2MS[15]=mWW2MS[20]*mWW2MS[15];
   mWW2MS[20]=mmt*mWW2MS[10];
   mWW2MS[21]= - 2./3.*mWW2MS[20] + 5*mWW2MS[6];
   mWW2MS[21]=mWW2MS[21]*mmt;
   mWW2MS[21]=4./3.*mWW2MS[9] + mWW2MS[21] - 6*mWW2MS[19] + 11./3.*
   mWW2MS[7];
   mWW2MS[21]=mmt*mWW2MS[21];
   mWW2MS[22]=4./3.*mWW2MS[11];
   mWW2MS[22]=mWW2MS[22]*pow(mmt,2);
   mWW2MS[17]=mWW2MS[22] - 20./3.*mWW2MS[17] + mWW2MS[21];
   mWW2MS[21]=mWW2MS[16] + 1;
   mWW2MS[16]=mWW2MS[16]*mWW2MS[21];
   mWW2MS[16]=mWW2MS[18] + mWW2MS[16];
   mWW2MS[16]=mWW2MS[4]*mWW2MS[17]*mWW2MS[16];
   mWW2MS[15]=2*mWW2MS[16] + mWW2MS[15];
   mWW2MS[15]=mWW2MS[4]*mWW2MS[15];
   mWW2MS[16]= - 1 - 2./3.*mWW2MS[6];
   mWW2MS[16]=mWW2MS[6]*mWW2MS[16];
   mWW2MS[16]=mWW2MS[20] + mWW2MS[16];
   mWW2MS[17]= - mWW2MS[9] + mWW2MS[7] - mWW2MS[19];
   mWW2MS[17]=mWW2MS[8]*mWW2MS[17];
   mWW2MS[17]=mWW2MS[17] - mWW2MS[11];
   mWW2MS[19]=mmZ*mWW2MS[10];
   mWW2MS[19]=8./3.*mWW2MS[19];
   mWW2MS[16]= - mWW2MS[19] - 13 + 4*mWW2MS[16] + 16./3.*mWW2MS[17];
   mWW2MS[16]=mWW2MS[16]*mWW2MS[18];
   mWW2MS[15]=mWW2MS[15] + mWW2MS[19] + mWW2MS[16];
   mWW2MS[15]=mWW2MS[1]*mWW2MS[15];
   mWW2MS[16]=mWW2MS[18]*mWW2MS[12];
   mWW2MS[17]=mWW2MS[12] - mWW2MS[16];
   mWW2MS[17]=mWW2MS[14]*mmZ*mWW2MS[17];
   mWW2MS[17]= - 31*mWW2MS[16] + 8*mWW2MS[17];
   mWW2MS[16]=mWW2MS[13]*mWW2MS[16];

      mWW2MSret = mWW2MS[15] - 4*mWW2MS[16] + 1./3.*mWW2MS[17];
      return mWW2MSret;
}
