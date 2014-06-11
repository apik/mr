#include <tt.hpp>
std::complex<long double> tt::mgl11(size_t nL, size_t nH)
{     
      
      
    std::complex<long double> mttgl[30];

    mttgl[1]=pow(SW,-1);
    mttgl[2]=pow(MMH,-1);
    mttgl[3]=pow(MMW,-1);
    mttgl[4]=pow(MMt,-1);
    mttgl[5]=Tsil::I2(MMH,MMt,MMt,mu2);
    mttgl[6]=Tsil::I2(0,0,MMt,mu2);
    mttgl[7]=Tsil::B(MMH,MMt,MMt,mu2);
    mttgl[8]=Tsil::A(MMH,mu2);
    mttgl[9]=Tsil::A(MMt,mu2);
    mttgl[10]=Tsil::B(0,0,MMt,mu2);
    mttgl[11]=Tsil::Beps(MMH,MMt,MMt,mu2);
    mttgl[12]=Tsil::Aeps(MMH,mu2);
    mttgl[13]=Tsil::Aeps(MMt,mu2);
    mttgl[14]=prot0ttHt->M(0);
    mttgl[15]=prottH0H->Vxzuv(0);
    mttgl[16]=prot0ttHt->Tuxv(0);
    mttgl[17]=protWt000->Tyzv(0);
    mttgl[18]=Mfin(MMH,MMt,MMt,0,MMt);
    mttgl[19]=Mfin(0,0,MMt,0,0);
    mttgl[20]=Mfin(0,0,0,MMt,0);
   mttgl[21]=mttgl[7] + 1 + mttgl[16];
   mttgl[22]=1./2.*mttgl[4];
   mttgl[21]=mttgl[22]*mttgl[21];
   mttgl[22]=mttgl[18] + mttgl[14];
   mttgl[21]=4*mttgl[15] + mttgl[21] - mttgl[22];
   mttgl[21]=MMH*mttgl[21];
   mttgl[23]=mttgl[4]*mttgl[9];
   mttgl[23]= - 17 - 7./2.*mttgl[23];
   mttgl[23]=mttgl[7]*mttgl[23];
   mttgl[21]=mttgl[23] + mttgl[21];
   mttgl[22]= - 8./3.*mttgl[15] + mttgl[22];
   mttgl[23]=2*MMt;
   mttgl[22]=mttgl[23]*mttgl[22];
   mttgl[23]=mttgl[5] - mttgl[8];
   mttgl[24]=1./6.*mttgl[4];
   mttgl[23]=mttgl[24]*mttgl[23];
   mttgl[24]=mttgl[17] - mttgl[11];
   mttgl[25]=mttgl[4]*pow(mttgl[9],2);
   mttgl[26]= - 5./3.*mttgl[13] - 1./3.*mttgl[9] + mttgl[25];
   mttgl[26]=mttgl[4]*mttgl[26];
   mttgl[21]= - 13./3. - 2*mttgl[16] + mttgl[23] + mttgl[22] + 
   mttgl[26] + 1./3.*mttgl[21] - 4./3.*mttgl[24];
   mttgl[21]=MMH*mttgl[21];
   mttgl[22]=mttgl[9]*mttgl[2];
   mttgl[23]=mttgl[8]*mttgl[2];
   mttgl[23]= - 4./3.*mttgl[23] + 7./2. - 13./3.*mttgl[22];
   mttgl[23]=mttgl[8]*mttgl[23];
   mttgl[24]=mttgl[22] + mttgl[17];
   mttgl[26]=mttgl[19] + mttgl[20];
   mttgl[26]= - 8./3.*mttgl[14] + 32*mttgl[2] - 1./3.*mttgl[26];
   mttgl[26]=MMt*mttgl[26];
   mttgl[24]=mttgl[26] + 52./3.*mttgl[7] - 8./3.*mttgl[11] + 13./3.*
   mttgl[16] - 17./6. + 4*mttgl[24];
   mttgl[24]=MMt*mttgl[24];
   mttgl[22]=23./2. - 36*mttgl[22];
   mttgl[22]=mttgl[9]*mttgl[22];
   mttgl[26]=7*mttgl[9] - mttgl[8];
   mttgl[26]=mttgl[7]*mttgl[26];
   mttgl[27]=mttgl[18]*pow(MMt,2);
   mttgl[28]=13./3.*mttgl[2];
   mttgl[28]=MMt*mttgl[28];
   mttgl[28]=3 + mttgl[28];
   mttgl[28]=mttgl[12]*mttgl[28];
   mttgl[29]=5./4.*mttgl[10] + 5./2.;
   mttgl[29]=MMt*mttgl[29];
   mttgl[29]=mttgl[29] + 4*mttgl[9];
   mttgl[29]=mttgl[10]*mttgl[29];
   mttgl[21]=1./3.*mttgl[29] + mttgl[28] + mttgl[21] - 8./3.*mttgl[27]
    + mttgl[24] + 2./3.*mttgl[26] - 13./6.*mttgl[6] - 11./3.*mttgl[5]
    + mttgl[23] + 65./6.*mttgl[13] + mttgl[22] - 2*mttgl[25];

      return mttgl[21]*mttgl[3]*pow(mttgl[1],2);
}
