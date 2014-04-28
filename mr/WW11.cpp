#include <WW.hpp>
std::complex<long double> WW::m11(size_t nG)
{     
      
      
    std::complex<long double> mWW[28];

    mWW[1]=pow(MMH,-1);
    mWW[2]=pow(MMW,-1);
    mWW[3]=double(nG);
    mWW[4]=prot00000->M(0);
    mWW[5]=Tsil::B(0,MMt,MMW,mu2);
    mWW[6]=pow(MMZ,-1);
    mWW[7]=Tsil::A(MMt,mu2);
    mWW[8]=pow(MMt,-1);
    mWW[9]=Tsil::Aeps(MMt,mu2);
    mWW[10]=prot00tt0->M(0);
    mWW[11]=prot00tt0->Tuxv(0);
    mWW[12]=1/(1 - pow(MMZ,-1)*MMW);
    mWW[13]=Tsil::B(0,0,MMW,mu2);
   mWW[14]= - 65./2. + 8*mWW[11];
   mWW[15]= - 1 + 1./3.*mWW[5];
   mWW[15]=mWW[5]*mWW[15];
   mWW[16]= - 4*mWW[9] - 5*mWW[7];
   mWW[17]=mWW[6]*mWW[16];
   mWW[14]=2./3.*mWW[17] + 1./3.*mWW[14] + 10*mWW[15];
   mWW[15]=mWW[6]*mWW[14];
   mWW[17]= - 4./3.*mWW[11] - 3*mWW[5];
   mWW[18]=mWW[6]*mWW[17];
   mWW[18]=32*mWW[1] + mWW[18];
   mWW[19]=mWW[6]*mWW[18];
   mWW[20]=MMt*mWW[10]*pow(mWW[6],2);
   mWW[19]=mWW[19] + 2./3.*mWW[20];
   mWW[19]=MMt*mWW[19];
   mWW[15]=2*mWW[19] - 4*mWW[10] + mWW[15];
   mWW[15]=MMt*mWW[15];
   mWW[19]= - 1 + mWW[5];
   mWW[20]= - 6*mWW[1] + mWW[8];
   mWW[20]=mWW[7]*mWW[20];
   mWW[19]=5./3.*mWW[19] + 2*mWW[20];
   mWW[19]=mWW[7]*mWW[19];
   mWW[20]=pow(mWW[7],2);
   mWW[21]=mWW[6]*mWW[20];
   mWW[19]=1./3.*mWW[21] + 2./3.*mWW[9] + mWW[19];
   mWW[21]=mWW[6]*mWW[19];
   mWW[22]=2./3. - mWW[13];
   mWW[23]=31./3. + 4*mWW[13];
   mWW[23]=mWW[3]*mWW[23];
   mWW[24]=mWW[3]*mWW[4];
   mWW[24]= - mWW[4] + mWW[24];
   mWW[24]=MMZ*mWW[24];
   mWW[25]=mWW[10]*MMZ;
   mWW[26]=mWW[9]*mWW[8];
   mWW[27]=1 + 2./3.*mWW[5];
   mWW[27]=mWW[5]*mWW[27];
   mWW[28]=mWW[5]*mWW[8];
   mWW[28]= - mWW[8] + mWW[28];
   mWW[28]=mWW[7]*mWW[28];
   mWW[15]=mWW[15] + 4*mWW[21] + 16./3.*mWW[28] + 4*mWW[27] + 16./3.*
   mWW[26] + 8./3.*mWW[25] + 16./3.*mWW[11] + 8./3.*mWW[24] + 4*mWW[22]
    + mWW[23];
   mWW[15]=mWW[12]*mWW[15];
   mWW[16]=mWW[2]*mWW[16];
   mWW[14]=mWW[14] + 2./3.*mWW[16];
   mWW[14]=mWW[2]*mWW[14];
   mWW[16]=mWW[2]*mWW[17];
   mWW[16]=mWW[18] + mWW[16];
   mWW[16]=mWW[2]*mWW[16];
   mWW[17]=mWW[6]*mWW[10];
   mWW[18]=mWW[2]*mWW[10];
   mWW[17]=mWW[17] + mWW[18];
   mWW[17]=MMt*mWW[2]*mWW[17];
   mWW[16]=mWW[16] + 2./3.*mWW[17];
   mWW[16]=MMt*mWW[16];
   mWW[14]=mWW[14] + 2*mWW[16];
   mWW[14]=MMt*mWW[14];
   mWW[16]= - mWW[3]*mWW[4];
   mWW[16]=mWW[4] + mWW[16];
   mWW[16]=MMZ*mWW[16];
   mWW[17]= - mWW[10]*MMZ;
   mWW[16]=mWW[16] + mWW[17];
   mWW[17]=mWW[2]*mWW[20];
   mWW[17]=mWW[19] + 1./3.*mWW[17];
   mWW[17]=mWW[2]*mWW[17];
   mWW[16]=2./3.*mWW[16] + mWW[17];

      return mWW[14] + mWW[15] + 4*mWW[16];
}
