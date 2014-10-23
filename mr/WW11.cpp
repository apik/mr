#include <WW.hpp>
std::complex<long double> WW::m11(size_t nL, size_t nH, int boson)
{     
      
      
    std::complex<long double> mWW[23], mWWret;

    mWW[1]=double(nH);
    mWW[2]=pow(CW,-1);
    mWW[3]=pow(MMH,-1);
    mWW[4]=pow(MMZ,-1);
    mWW[5]=pow(SW,-1);
    mWW[6]=Tsil::B(0,MMt,MMW,mu2);
    mWW[7]=Tsil::A(MMt,mu2);
    mWW[8]=pow(MMt,-1);
    mWW[9]=Tsil::Aeps(MMt,mu2);
    mWW[10]=prot00tt0->M(0);
    mWW[11]=prot00tt0->Tuxv(0);
    mWW[12]=double(nL);
    mWW[13]=std::real(Tsil::B(0,0,MMW,mu2));
    mWW[14]=prot00000->M(0);
   mWW[15]= - 1 + 1./3.*mWW[6];
   mWW[16]=10*mWW[6];
   mWW[15]=mWW[15]*mWW[16];
   mWW[16]= - 65./2. + 8*mWW[11];
   mWW[15]=mWW[15] + 1./3.*mWW[16];
   mWW[15]=mWW[15]*MMt;
   mWW[16]=mWW[6] - 1;
   mWW[17]=5*mWW[7];
   mWW[18]=mWW[16]*mWW[17];
   mWW[18]=mWW[18] + 2*mWW[9];
   mWW[19]=pow(mWW[7],2);
   mWW[20]=mWW[19]*mWW[8];
   mWW[18]=2*mWW[20] + 1./3.*mWW[18];
   mWW[20]=pow(MMt,2);
   mWW[20]= - 4*mWW[20] + 3*mWW[19];
   mWW[21]=16*mWW[3];
   mWW[20]=mWW[20]*mWW[21];
   mWW[15]= - mWW[20] + mWW[15] + 4*mWW[18];
   mWW[18]=pow(mWW[2],2);
   mWW[20]=pow(mWW[5],2);
   mWW[21]=mWW[20] + mWW[18];
   mWW[15]=mWW[15]*mWW[21];
   mWW[21]= - mWW[18] - 1;
   mWW[18]=mWW[18]*mWW[21];
   mWW[18]= - mWW[20] + mWW[18];
   mWW[21]=MMt*mWW[10];
   mWW[22]= - 2./3.*mWW[21] + 3*mWW[6] + 4./3.*mWW[11];
   mWW[22]=mWW[22]*MMt;
   mWW[17]=mWW[17] + 4*mWW[9];
   mWW[17]=mWW[22] + 1./3.*mWW[17];
   mWW[17]=mWW[17]*MMt;
   mWW[17]=mWW[17] - 2./3.*mWW[19];
   mWW[17]=mWW[4]*mWW[17]*mWW[18];
   mWW[15]=2*mWW[17] + mWW[15];
   mWW[15]=mWW[4]*mWW[15];
   mWW[17]=1 + 2./3.*mWW[6];
   mWW[17]=mWW[6]*mWW[17];
   mWW[17]=mWW[21] - mWW[17];
   mWW[16]=mWW[7]*mWW[16];
   mWW[16]=mWW[9] + mWW[16];
   mWW[16]=mWW[8]*mWW[16];
   mWW[16]=mWW[16] + mWW[11];
   mWW[18]=MMZ*mWW[10];
   mWW[18]=8./3.*mWW[18];
   mWW[16]=mWW[18] + 13 - 4*mWW[17] + 16./3.*mWW[16];
   mWW[16]=mWW[16]*mWW[20];
   mWW[15]=mWW[15] - mWW[18] + mWW[16];
   mWW[15]=mWW[1]*mWW[15];
   mWW[16]=mWW[14]*MMZ*mWW[12];
   mWW[16]=8*mWW[16];
   mWW[17]=31*mWW[12] + mWW[16];
   mWW[17]=mWW[17]*mWW[20];
   mWW[16]= - mWW[16] + mWW[17];
   mWW[17]=mWW[13]*mWW[12]*mWW[20];

      mWWret = mWW[15] + 1./3.*mWW[16] + 4*mWW[17];
      return mWWret;
}
